""" This file is part of SRC.

    SRC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SRC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SRC.  If not, see <https://www.gnu.org/licenses/>.

    Copyright 2020 David T. Mildebrath, Logan A. Smith
"""
import networkx as nx
import gurobipy as grb
import prepr
import clique
import heur
import time
import math


def src_BnC(
    G,
    max_time=3600,
    clique_method="ostergard",
    fix_method="clique",
    heuristic_path_fixing=False,
    gurobi_verbose=False,
    verbose=False,
    add_clique_cuts=False,
    print_computed_solution=False,
    break_symmetry=True,
    heur_bound=None,
    num_heur_trials=None,
):
    """ Compute src(G) for the specified graph using branch-and-cut.

    Args:
        G (networkx.Graph):
            A simple, connected graph.
        max_time (int, optional):
            Maximum computational time allowed in seconds. Default is 3600.
        clique_method (str, optional):
            Method used to find cliques in the auxiliary graph H. Must be one
            of the following:
            - `brute`: Enumerate all cliques by brute force and select the
              largest one.
            - `ostergard`: Find the largest clique using the branch-and-bound
              method of Ostergard (aka cliquer).
            - `ILS_VND`: Heuristically look for the largest clique using the
              ILS_VND heuristic.
            Default is `ostergard`.
        fix_method (str, optional):
            Method used to fix an initial set of edges in the model. Must be
            one of the following:
            - `lusp`: Fix every edge in a longest unique shortest path.
            - `clique`: Fix every edge corresponding to a clique in the
              auxiliary graph H.
            - `none`: No initial set of edges is fixed.
            Default is `clique`. 
        heuristic_path_fixing (bool, optional):
            If True, select one shortest path between each pair of vertices,
            and only require that path to be rainbow colored in the result.
            This produces a valid strong rainbow coloring of G, which may not
            be minimum. Default is False.
        gurobi_verbose (bool, optional):
            If True, the Gurobi output flag is set to True. Default is False.
        verbose (bool, optional):
            If True, displays verbose output. Default is False.
        add_clique_cuts (bool, optional):
            If True, the Gurobi LazyConstraint parameter is set to True.
            Default is False.
        print_computed_solution (bool, optional):
            If True, prints the solution. Default is False.
        break_symmetry (bool, optional):
            If True, add symmetry breaking constraints of the form z_k >=
            z_{k+1} for all colors k. Default is True.
        heur_bound (int, optional):
            Valid upper bound on src(G). If none is specified, an upper bound
            is obtained with the heuristic. Default is None (no bound).
        num_heur_trials (int, optional):
            Number of heuristic trials to run. Ignored if heur_bound is
            specified. Default is ceil(n / 5).


    Returns:
        float or None:
            Returns src(G) is computed successfully, otherwise None.
    """
    start_time = time.time()
    prepr.presolve(G)

    # Build C/C++ graphs as necessary
    if clique_method == "ILS_VND":
        G.Hcpp = clique.build_cppgraph(G)
        G.Hcpp.sort() # Is this necessary?
    elif clique_method == "ostergard":
        G.Hc = clique.build_c_graph(G)

    # Compute an initial set of edges to fix
    if fix_method.lower() != "none":
        G.init_fix = fix_initial_edge_set(
            G,
            fix_method=fix_method,
            clique_method=clique_method,
        )

    # If desired, fix a random path for each pair of vertices
    # (used for producing a heuristic solution)
    init_paths = None
    if heuristic_path_fixing:
        init_paths = {
            (u, v): (nx.shortest_path(G, source=u, target=v),)
            for (u, v) in G.vertex_pairs
        }

    # Run the heuristic
    if heur_bound is None:
        if num_heur_trials is None:
            num_heur_trials = math.ceil(G.order() / 5)

        start = time.time()
        if verbose:
            print("Beginning heuristic...", end="", flush=True)
        _, bound = heur.path_fixing_heuristic(G, num_trials=num_heur_trials)
        if verbose:
            print("done.")
        end = time.time()
        K = bound
        if verbose:
            print(f"Heuristic provides an upper bound of {bound}, in {end - start:.4f} seconds.")
    else:
        K = heur_bound 

    # Build the model
    model = build_model(
        G,
        K,
        break_symmetry=break_symmetry,
        init_paths=init_paths,
        gurobi_verbose=gurobi_verbose,
        add_clique_cuts=add_clique_cuts,
        gurobi_break_symmetry=False,
    )

    # Fix the initial set of edges
    if fix_method.lower() != "none":
        fix_init_vars(G, model)
        if verbose:
            print(f"Number of edges fixed by init_fix: {len(G.init_fix)}")

    if add_clique_cuts:
        # This is used in the cut callback functions
        model._clique_method = clique_method

    if gurobi_verbose:
        print("\n################## BEGIN GUROBI OUTPUT #########################")

    time_elapsed = time.time() - start_time
    time_remaining = max(0, max_time - time_elapsed)
    model.Params.TimeLimit = time_remaining

    if add_clique_cuts:
        model.optimize(clique_callback)
    else:
        model.optimize()

    if gurobi_verbose:
        print("################### END GUROBI OUTPUT ##########################\n")

    ## THIS IS IMPORTANT! (Leaks otherwise)
    if clique_method == "ostergard":
        clique.free_graph(G.Hc)

    if model.Status == grb.GRB.Status.OPTIMAL:
        if print_computed_solution:
            print_solution(G, model)
        return model.objval


def src_bottom_up(
    G,
    max_time=3600,
    clique_method="ostergard",
    fix_method="clique",
    heuristic_path_fixing=False,
    gurobi_verbose=False,
    verbose=False,
    add_clique_cuts=False,
    print_computed_solution=False,
    break_symmetry=True,
):
    """ Compute src(G) for the specified graph using the "bottom-up" approach.

    Args:
        G (networkx.Graph):
            A simple, connected graph.
        max_time (int, optional):
            Maximum computational time allowed in seconds. Default is 3600.
        clique_method (str, optional):
            Method used to find cliques in the auxiliary graph H. Must be one
            of the following:
            - `brute`: Enumerate all cliques by brute force and select the
              largest one.
            - `ostergard`: Find the largest clique using the branch-and-bound
              method of Ostergard (aka cliquer).
            - `ILS_VND`: Heuristically look for the largest clique using the
              ILS_VND heuristic.
            Default is `ostergard`.
        fix_method (str, optional):
            Method used to fix an initial set of edges in the model. Must be
            one of the following:
            - `lusp`: Fix every edge in a longest unique shortest path.
            - `clique`: Fix every edge corresponding to a clique in the
              auxiliary graph H.
            - `none`: No initial set of edges is fixed.
            Default is `clique`. 
        heuristic_path_fixing (bool, optional):
            If True, select one shortest path between each pair of vertices,
            and only require that path to be rainbow colored in the result.
            This produces a valid strong rainbow coloring of G, which may not
            be minimum. Default is False.
        gurobi_verbose (bool, optional):
            If True, the Gurobi output flag is set to True. Default is False.
        verbose (bool, optional):
            If True, displays verbose output. Default is False.
        add_clique_cuts (bool, optional):
            If True, the Gurobi LazyConstraint parameter is set to True.
            Default is False.
        print_computed_solution (bool, optional):
            If True, prints the solution. Default is False.
        break_symmetry (bool, optional):
            If True, add symmetry breaking constraints of the form z_k >=
            z_{k+1} for all colors k. Default is True.

    Returns:
        float or None:
            Returns src(G) is computed successfully, otherwise None.
    """
    start_time = time.time()
    prepr.presolve(G)

    # Build C/C++ graphs as necessary
    if clique_method == "ILS_VND":
        G.Hcpp = clique.build_cppgraph(G)
        G.Hcpp.sort() # Is this necessary?
    elif clique_method == "ostergard":
        G.Hc = clique.build_c_graph(G)

    # Compute an initial set of edges to fix
    if fix_method.lower() != "none":
        G.init_fix = fix_initial_edge_set(
            G,
            fix_method=fix_method,
            clique_method=clique_method,
        )

    # If desired, fix a random path for each pair of vertices
    # (used for producing a heuristic solution)
    init_paths = None
    if heuristic_path_fixing:
        init_paths = {
            (u, v): (nx.shortest_path(G, source=u, target=v),)
            for (u, v) in G.vertex_pairs
        }

    K = max(len(G.init_fix), G.diameter)

    model = build_model(
        G,
        K,
        break_symmetry=break_symmetry,
        init_paths=init_paths,
        gurobi_verbose=gurobi_verbose,
        add_clique_cuts=add_clique_cuts,
        gurobi_break_symmetry=False,
    )

    # Fix the initial set of edges
    if fix_method.lower() != "none":
        fix_init_vars(G, model)
        if verbose:
            print(f"Number of edges fixed by init_fix: {len(G.init_fix)}")

    if add_clique_cuts:
        # This is used in the cut callback functions
        model._clique_method = clique_method

    i = 0
    sol_found = False
    while not sol_found and time.time() - start_time < max_time:
        time_elapsed = time.time() - start_time
        model.Params.TimeLimit = max(max_time - time_elapsed, 0)
        if gurobi_verbose:
            print(f"\n\n################################################################")
            print(f"### Beginning Iteration: {i} #####################################")
            print(f"################################################################")
            print("\n################## BEGIN GUROBI OUTPUT #########################")
        elif verbose:
            print("Beginning iteration", i)

        if add_clique_cuts:
            model.optimize(bottom_up_clique_callback)
        else:
            model.optimize(bottom_up_no_cut_callback)

        if gurobi_verbose:
            print("################### END GUROBI OUTPUT ##########################\n")

        sol_found = model.SolCount > 0
        if not sol_found:
            i += 1
            K += 1
            add_color(model, break_symmetry=break_symmetry)

    ## THIS IS IMPORTANT! (Leaks otherwise)
    if clique_method == "ostergard":
        clique.free_graph(G.Hc)

    if sol_found:
        if print_computed_solution:
            print_solution(G, model)
        return model.objval


def build_model(
    G,
    K,
    break_symmetry=True,
    init_paths=None,
    gurobi_verbose=False,
    add_clique_cuts=False,
    gurobi_break_symmetry=False,
    timeout=None,
):
    """ Build the IP model for computing src(G).

    Args:
        G (networkx.Graph):
            A simple, connected graph. Must have attributes for `ordered_edges`
            and `vertex_pairs`.
        K (int):
            The number of colors in the palette.
        break_symmetry (bool, optional):
            If True, add symmetry breaking constraints of the form z_k >=
            z_{k+1} for all colors k. Default is True.
        init_paths (dict, optional):
            Dictionary mapping vertex pairs (u, v) to list of shortest (u, v)
            paths in G. If init_paths[u, v] exists, then only those paths will
            be added to the resulting IP model. If it does not, then _all_
            shortest (u, v) paths in G will be included. Default is an empty
            dictionary (i.e. all shortest (u, v) paths are included in the
            model).
        gurobi_verbose (bool, optional):
            If True, the Gurobi output flag is set to True. Default is False.
        add_clique_cuts (bool, optional):
            If True, the Gurobi LazyConstraint parameter is set to True.
            Default is False.
        gurobi_break_symmetry (bool, optional):
            If True, sets the Gurobi Symmetry parameter to 2. Default is False.
        timeout (int, optional):
            If specified, sets the Gurobi TimeLimit parameter (units of
            seconds). Default is None (no time limit).

    Returns:
        gurobipy.Model:
            The resulting IP model.

    Notes:
        - This function manually sets the LP solution method to be dual
          simplex, in order to prevent the barrier method from being used. In
          practice the dual simplex method is much more effective for this
          problem.
    """
    required_attrs = ["vertex_pairs", "ordered_edges"]
    for attr in required_attrs:
        if not hasattr(G, attr):
            raise RuntimeError(
                f"Graph must have a '{attr}' attribute before building model"
            )

    if init_paths is None:
        init_paths = dict()

    model = grb.Model()

    K = range(K)
    x = model.addVars(G.ordered_edges, K, vtype=grb.GRB.BINARY, name="x")
    z = model.addVars(K, vtype=grb.GRB.BINARY, name="z")

    obj = grb.quicksum(z[k] for k in K)
    model.setObjective(obj, sense=grb.GRB.MINIMIZE)

    """ Add the basic logical constraints """
    model._sum_one = grb.tupledict()
    for (i, j) in G.ordered_edges:
        model._sum_one[i, j] = model.addConstr(grb.quicksum(x[i, j, k] for k in K) == 1)
        for k in K:
            model.addConstr(x[i, j, k] <= z[k])

    """ Now the tricky constraints """
    y = grb.tupledict()

    if break_symmetry:
        for k in range(len(K) - 1):
            model.addConstr(z[k] >= z[k + 1])

    for (u, v) in G.vertex_pairs:
        if (u, v) in init_paths:
            paths = init_paths[u, v]
        else:
            paths = nx.all_shortest_paths(G, source=u, target=v)
            paths = list(paths)
        for (i, p) in enumerate(paths):
            y[u, v, i] = model.addVar(vtype=grb.GRB.BINARY, name=f"y[{u},{v},{i}]")
            edges_in_path = prepr.path_to_edge_tuples(p, G.ordered_edges)
            for k in K:
                model.addConstr(
                    grb.quicksum(x[n1, n2, k] for (n1, n2) in edges_in_path)
                    + (len(p) - 1) * y[u, v, i] <= len(p)
                )
        if len(paths) == 1:
            model.addConstr(y[u, v, i] == 1)
        else:
            model.addConstr(grb.quicksum(y[u, v, i] for i in range(len(paths))) >= 1)

    model._K = K
    model._x = x
    model._y = y
    model._z = z
    model._G = G

    model.Params.OutputFlag = gurobi_verbose

    if add_clique_cuts:
        model.Params.LazyConstraints = 1

    if gurobi_break_symmetry:
        model.Params.Symmetry = 2

    if timeout is not None:
        model.Params.TimeLimit = timeout

    ## !!!
    model.Params.Method = 1 # Dual simplex--no barrier!

    model.update()
    return model


def add_color(model, break_symmetry=True):
    """ Add a color to the color palette (for bottom up).

    Args:
        model (gurobipy.Model):
            IP model. Must have attached graph, variables, etc.
        break_symmetry (bool, optional):
            If True, adds a symmetry breaking constraint of the form z_k <=
            z_{k+1} for the new colors.
    """
    # Add color to palette
    model._K = range(len(model._K) + 1)
    K = len(model._K) # Number of colors

    # Add new z variable and update objective
    model._z[K-1] = model.addVar(vtype=grb.GRB.BINARY, obj=1, name=f"z[{K-1}]")

    # Add new x variables
    for (u, v) in model._G.ordered_edges:
        model._x[u, v, K-1] = model.addVar(
            vtype=grb.GRB.BINARY, name=f"x[{u},{v},{K-1}]"
        )

    # Add new variable upper bound constraint
    for (u, v) in model._G.ordered_edges:
        model.addConstr(model._x[u, v, K-1] <= model._z[K-1])

    # Add new symmetry breaking constraint
    if break_symmetry:
        model.addConstr(model._z[K-2] >= model._z[K-1])

    # Add |P| new "big-M" path constraints
    for (u, v) in model._G.vertex_pairs:
        paths = list(nx.all_shortest_paths(model._G, source=u, target=v))
        for (i, p) in enumerate(paths):
            edges_in_path = prepr.path_to_edge_tuples(p, model._G.ordered_edges)
            model.addConstr(
                grb.quicksum(model._x[n1, n2, K-1] for (n1, n2) in edges_in_path)
                + (len(p) - 1) * model._y[u, v, i] <= len(p)
            )

    # Add the new x variables to the sum-one constraints
    for (u, v) in model._G.ordered_edges:
        model.chgCoeff(model._sum_one[u, v], model._x[u, v, K-1], 1)

    # Keep the init_fix up to date
    if hasattr(model._G, "init_fix") and K > len(model._G.init_fix):
        for (u, v) in model._G.init_fix:
            model.addConstr(model._x[u, v, K-1] == 0)

    model.update()


def clique_callback(model, where):
    if where != grb.GRB.Callback.MIPNODE :
        return
    if model.cbGet(grb.GRB.Callback.MIPNODE_STATUS) != grb.GRB.OPTIMAL:
        return
    # At this point, we have an optimal LP solution at this node
    G = model._G

    for k in model._K:
        vbs = {(u, v): model._x[u, v, k] for (u, v) in G.ordered_edges}
        xk = model.cbGetNodeRel(vbs)
        totalWeight = 0
        if sum(xk.values()) <= 1:
            continue
        if model._clique_method == "brute":
            totalWeight, constr = clique.enumerate_best_clique(G.aux_cut_graph, xk)
        elif model._clique_method == "ostergard":
            totalWeight, constr = clique.max_clique_ostergard(G, weights=xk, scale=10000)
        elif model._clique_method == "ILS_VND":
            totalWeight, constr = clique.max_clique_ILS_VND(G, weights=xk, scale=10000, iterations=10)
        if totalWeight > 1 + 1e-8:
            model.cbLazy(grb.quicksum(model._x[u,v,k] for (u,v) in constr) <= 1)

    return


def bottom_up_clique_callback(model, where):
    if where == grb.GRB.Callback.MIPSOL :
        model.terminate()
    if where != grb.GRB.Callback.MIPNODE :
        return
    if model.cbGet(grb.GRB.Callback.MIPNODE_STATUS) != grb.GRB.OPTIMAL:
        return
    # At this point, we have an optimal LP solution at this node
    G = model._G

    for k in model._K :
        vbs = {(u, v): model._x[u, v, k] for (u, v) in G.ordered_edges}
        xk = model.cbGetNodeRel(vbs)
        totalWeight = 0
        if sum(xk.values()) <= 1 :
            continue
        if model._clique_method == "brute" :
            totalWeight, constr = clique.enumerate_best_clique(G.aux_cut_graph, xk)
        elif model._clique_method == "ostergard" :
            totalWeight, constr = clique.max_clique_ostergard(G, weights=xk, scale=10000)
        elif model._clique_method == "ILS_VND" :
            totalWeight, constr = clique.max_clique_ILS_VND(G, weights=xk, scale=10000, iterations=10)
        if totalWeight > 1 + 1e-8 :
            model.cbLazy(grb.quicksum(model._x[u,v,k] for (u,v) in constr) <= 1)
    return


def bottom_up_no_cut_callback(model, where):
    if where == grb.GRB.Callback.MIPSOL :
        model.terminate()
    return


def fix_initial_edge_set(
    G,
    fix_method="clique",
    clique_method="ostergard",
    ils_vnd_iterations=10,
    clique_scale=10000,
):
    """ Compute a set of edges which must all be different colors in any valid
    strong rainbow coloring of G.

    Args:
        G (Graph):
            The graph.
        fix_method (str, optional):
            The method used to compute the edge set. Must be one of the
            following:
            - `lusp`: Fix every edge in a longest unique shortest path.
            - `clique`: Fix every edge corresponding to a clique in the
              auxiliary graph H.
            Default is `clique`. 
        clique_method (str, optional):
            Method used to find cliques in the auxiliary graph H. Ignored
            unless fix_method is `clique`. Must be one of the following:
            - `brute`: Enumerate all cliques by brute force and select the
              largest one.
            - `ostergard`: Find the largest clique using the branch-and-bound
              method of Ostergard (aka cliquer).
            - `ILS_VND`: Heuristically look for the largest clique using the
              ILS_VND heuristic.
            Default is `ostergard`.
        ils_vnd_iterations (int, optional):
            Number of ILS_VND iterations to perform. Ignored unless
            clique_method is `ILS_VND`. Default is 10.
        clique_scale (int, optional):
            Amount to scale edge weight by before casting to integers in the
            ostergard/ILS_VND methods (ignored if neither of these methods are
            used). Default is 10000.

    Returns:
        list:
            List of edges which must be different colors in any valid strong
            rainbow coloring.
    """
    if fix_method not in ("lusp", "clique", "none"):
        raise ValueError(f"Unrecognized value `{init_fix}` passed for fix_method")
    if fix_method == "clique":
        if clique_method not in ("brute", "ostergard", "ILS_VND"):
            raise ValueError(f"Unrecognized value `{clique_method}` passed for clique_method")

    if fix_method == "lusp":
        initial_edges = prepr.find_longest_unique_shortest_path(G)
    elif fix_method == "clique":
        if clique_method == "brute":
            _, initial_edges = clique.enumerate_best_clique(G.aux_cut_graph)
        elif clique_method == "ostergard" :
            _, initial_edges = clique.max_clique_ostergard(
                G,
                weights=None,
                scale=clique_scale,
            )
        elif clique_method == "ILS_VND" :
            _, initial_edges = clique.max_clique_ILS_VND(
                G,
                weights=None,
                scale=clique_scale,
                iterations=ils_vnd_iterations,
            )
    return initial_edges


def graph_statistics(G, verbose=True):
    prepr.presolve(G)

    num_total_paths = sum(G.path_counts.values())
    num_remaining_paths = sum(G.path_counts[u, v] for (u, v) in G.vertex_pairs)

    G.H_density = 2 * G.size() / (G.order() * (G.order() - 1))

    # Compute clique bound
    G.Hc = clique.build_c_graph(G)
    _, G.init_fix = clique.max_clique_ostergard(G, weights=None, scale=10000)
    clique.free_graph(G.Hc)

    if verbose:
        """ Print some graph diagnostics """
        uniques = sum(1 for (u, v) in G.vertex_pairs if G.path_counts[u, v] == 1)
        num_pairs = int(G.order() * (G.order() - 1) / 2)

        print(f"n: {G.order()}")
        print(f"m: {G.size()}")
        print(f"Graph diameter: {G.diameter}")
        print(f"Graph w' bound: {len(G.init_fix)}")
        print(f"Edge density of auxiliary graph: {100*G.H_density:.2f}%")
        print(f"Paths in G**: {num_total_paths}")
        print(f"Paths in G** (after elimination): {num_remaining_paths} ({100*num_remaining_paths/num_total_paths:.2f}%)")
        print(f"Non-adjacent vertex pairs in G: {num_pairs}")
        print(f"Non-adjacent vertex pairs in G (after elimination): {len(G.vertex_pairs)}")
        print(f"Non-adjacent vertex pairs in G (after elimination) with unique shortest paths: {uniques}")
        print(f"Mean paths/vertex pair (all paths): {num_total_paths/num_pairs:.4f}")
        print(f"Mean paths/vertex pair (after elimination): {num_remaining_paths/len(G.vertex_pairs):.4f}")
        print()
        print("** Does not count single edges as paths")


def fix_init_vars(G, model):
    if hasattr(G, "init_fix"):
        for (k, (u, v)) in enumerate(G.init_fix):
            for kp in model._K:
                if kp == k:
                    model.addConstr(model._x[u, v, kp] == 1)
                else:
                    model.addConstr(model._x[u, v, kp] == 0)


def print_solution(G, model):
    colors = dict()
    for (u, v, k) in model._x.keys():
        if model._x[u, v, k].X > 0.5:
            if k in colors:
                colors[k].append((u, v))
            else:
                colors[k] = [(u, v)]
    for (i, k) in enumerate(list(colors.keys())):
        print(f"Color {i+1}:", ", ".join(map(str, colors[k])))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Compute the strong rainbow connection number a graph"
    )
    parser.add_argument("filename", help="Path to DIMACS file")
    parser.add_argument(
        "method", type=str, choices=["bu", "bnc"],
        help="Method to compute the src. Must be either `bnc` for "
        "branch-and-cut or `bu` for bottom up."
    )
    parser.add_argument(
        "-g", "--grb-verbose", action="store_true", default=False,
        help="Display Gurobi output. Default is False."
    )
    parser.add_argument(
        "-k", "--colors", type=int,
        help="Number of colors to use in the palette (ignored if using bottom up)."
    )
    parser.add_argument(
        "-b", "--no-symmetry", action="store_true",
        help="Do not add color symmetry breaking constraints."
    )
    parser.add_argument(
        "-x", "--init-fix", type=str, default="clique",
        choices=["clique", "lusp", "none"],
        help="Type of initial edge fixing to use. Default is `clique`."
    )
    parser.add_argument(
        "-c", "--clique-cuts", action="store_true",
        help="Use clique cuts."
    )
    parser.add_argument(
        "-s", "--no-skip-pairs", action="store_true",
        help="Do not skip dominated pairs of paths."
    )
    parser.add_argument(
        "-z", "--clique-method", type=str, default="ostergard",
        choices=["ostergard", "brute", "none"],
        help="Method to use for computing cliques in auxiliary graph H. Default is `ostergard`"
    )
    parser.add_argument(
        "-t", "--timeout", type=float, default=3600,
        help="Timeout for the given method, in seconds. Default is 3600s."
    )
    parser.add_argument(
        "-p", "--print-solution", action="store_true",
        help="Display the resulting edge coloring"
    )
    parser.add_argument(
        "-n", "--heur-trials", type=int,
        help="Number of heuristic trials to run (only used when method=bnc). "
        "Default is the number of vertices in the graph divided by 5."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Turn on output from various methods",
    )
    args = parser.parse_args()

    G = prepr.build_graph_from_file(args.filename)

    t0 = time.time()
    if args.method == "bu":
        srcg = src_bottom_up(
            G,
            max_time=args.timeout,
            clique_method=args.clique_method,
            fix_method=args.init_fix,
            heuristic_path_fixing=False,
            gurobi_verbose=args.grb_verbose,
            verbose=args.verbose,
            add_clique_cuts=args.clique_cuts,
            print_computed_solution=args.print_solution,
            break_symmetry=not args.no_symmetry,
        )
    else:
        srcg = src_BnC(
            G,
            max_time=args.timeout,
            clique_method=args.clique_method,
            fix_method=args.init_fix,
            heuristic_path_fixing=False,
            gurobi_verbose=args.grb_verbose,
            verbose=args.verbose,
            add_clique_cuts=args.clique_cuts,
            print_computed_solution=args.print_solution,
            break_symmetry=not args.no_symmetry,
            heur_bound=args.colors,
            num_heur_trials=args.heur_trials,
        )
    if srcg is not None:
        print(f"src = {srcg:.1f}")
    else:
        print("Could not compute src")
    tf = time.time()
    print(f"Total time: {tf-t0:.2f}s")
