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


def src_BnC(G, settings):
    """ Compute src(G) for the specified graph using branch-and-cut.

    Args:
        G (networkx.Graph):
            A simple, connected graph.
        settings (dict):
            Dictionary of settings to control the solution. See the
            documentation for `src.check_settings` for details.

    Returns:
        float or None:
            Returns src(G) is computed successfully, otherwise None.

    Notes:
        - This function requires the settings dict to contain a key `timeout`
          pointing to a float equal to the maximum run time in seconds.
    """
    start_time = time.time()
    max_time = settings["timeout"]
    check_settings(G, settings)
    graph_preprocessing(G, settings)

    # Run the heuristic
    heur_bound = None
    if "heur_bound" in settings:
        heur_bound = settings["heur_bound"]

    if heur_bound is None:
        num_heur_trials = math.ceil(G.order() / 5)
        if "num_heur_trials" in settings and settings["num_heur_trials"] is not None:
            num_heur_trials = settings["num_heur_trials"]
        start = time.time()
        _, bound = heur.path_fixing_heuristic(G, num_trials=num_heur_trials)
        end = time.time()
        settings["K"] = bound
        if settings["verbose"]:
            print(f"Heuristic provides an upper bound of {bound}, in {end - start:.4f} seconds.")
    else:
        settings["K"] = heur_bound

    time_elapsed = time.time() - start_time
    settings["timeout"] = max(0, max_time - time_elapsed)
    model = build_model(G, settings)
    model_preprocessing(G, model, settings)

    if settings["GRB_verbose"]:
        print("\n################## BEGIN GUROBI OUTPUT #########################")

    if settings["add_clique_cuts"] :
        model.optimize(clique_callback)
    else :
        model.optimize()

    if settings["GRB_verbose"]:
        print("################### END GUROBI OUTPUT ##########################\n")

    if settings["print_solution"]:
        print_solution(G, model)

    ## THIS IS IMPORTANT! (Leaks otherwise)
    if settings["clique_method"] == "ostergard":
        clique.free_graph(G.Hc)

    if model.Status == grb.GRB.Status.OPTIMAL :
        return model.objval


def src_bottom_up(G, settings):
    """ Compute src(G) for the specified graph using the "bottom-up" approach.

    Args:
        G (networkx.Graph):
            A simple, connected graph.
        settings (dict):
            Dictionary of settings to control the solution. See the
            documentation for `src.check_settings` for details.

    Returns:
        float or None:
            Returns src(G) is computed successfully, otherwise None.

    Notes:
        - This function requires the settings dict to contain a key `timeout`
          pointing to a float equal to the maximum run time in seconds.
    """
    start_time = time.time()
    max_time = settings["timeout"]
    check_settings(G, settings)
    graph_preprocessing(G, settings)
    settings["K"] = max(len(G.init_fix), G.diameter)

    model = build_model(G, settings)
    model_preprocessing(G, model, settings)

    i = 0
    sol_found = False
    while not sol_found and time.time() - start_time < max_time:
        time_elapsed = time.time() - start_time
        model.Params.TimeLimit = max(max_time - time_elapsed, 0)
        if settings["GRB_verbose"]:
            print(f"\n\n################################################################")
            print(f"### Beginning Iteration: {i} #####################################")
            print(f"################################################################")
            print("\n################## BEGIN GUROBI OUTPUT #########################")
        elif settings["verbose"]:
            print("Beginning iteration", i)

        if settings["add_clique_cuts"] :
            model.optimize(bottom_up_clique_callback)
        else :
            model.optimize(bottom_up_no_cut_callback)

        if settings["GRB_verbose"]:
            print("################### END GUROBI OUTPUT ##########################\n")

        if model.SolCount == 0 : # No solution found yet
            i += 1
            settings["K"] += 1
            add_color(model) # No need to call model_preprocesing again

        else :
            sol_found = True
            if settings["print_solution"]:
                print_solution(G, model)

    ## THIS IS IMPORTANT! (Leaks otherwise)
    if settings["clique_method"] == "ostergard":
        clique.free_graph(G.Hc)

    if sol_found :
        return model.objval


def build_model(G, settings):
    """ Build the IP model for computing src(G).

    Args:
        G (networkx.Graph):
            A simple, connected graph.
        settings (dict):
            Dictionary of settings to control the solution. See the
            documentation for `src.check_settings` for details.

    Returns:
        gurobipy.Model:
            The resulting IP model.

    Notes:
        - This function manually sets the LP solution method to be dual
          simplex, in order to prevent the barrier method from being used. In
          practice the dual simplex method is much more effective for this
          problem.
        - The user may specifiy an `init_paths` key in the settings dict. This
          key should point to another dict which maps (u,v) vertex pairs to
          list of paths which the model may select from. If `init_paths` is not
          specified (or is missing a vertex pair), the set of all shortest
          (u,v) is explicitly computed using networkx.
    """
    required_attrs = ["vertex_pairs", "ordered_edges"]
    for attr in required_attrs:
        if not hasattr(G, attr):
            raise RuntimeError(
                f"Graph must have a '{attr}' attribute before building model"
            )

    model = grb.Model()

    K = range(settings["K"])
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

    break_symmetry = settings["add_symmetry_breaking"]
    if break_symmetry:
        for k in range(len(K) - 1):
            model.addConstr(z[k] >= z[k + 1])

    init_paths = settings["init_paths"]
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
    model._settings = settings

    if not settings["GRB_verbose"]:
        model.Params.OutputFlag = settings["GRB_verbose"]

    if settings["add_clique_cuts"] :
        model.Params.LazyConstraints = 1

    if settings["GRB_symmetry_breaking"] :
        model.Params.Symmetry = 2

    if settings["timeout"] is not None :
        model.Params.TimeLimit = settings["timeout"]

    ## !!!
    model.Params.Method = 1 # Dual simplex--no barrier!

    model.update()
    return model


def add_color(model):
    """ Add a color to the color palette (for bottom up)
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

    for k in model._K :
        vbs = {(u, v): model._x[u, v, k] for (u, v) in G.ordered_edges}
        xk = model.cbGetNodeRel(vbs)
        totalWeight = 0
        if sum(xk.values()) <= 1 :
            continue
        if model._settings["clique_method"] == "brute" :
            totalWeight, constr = clique.enumerate_best_clique(G.aux_cut_graph, xk)
        elif model._settings["clique_method"] == "ostergard" :
            totalWeight, constr = clique.max_clique_ostergard(G, weights=xk, scale=10000)
        elif model._settings["clique_method"] == "ILS_VND" :
            totalWeight, constr = clique.max_clique_ILS_VND(G, weights=xk, scale=10000, iterations=10)
        if totalWeight > 1 + 1e-8 :
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
        if model._settings["clique_method"] == "brute" :
            totalWeight, constr = clique.enumerate_best_clique(G.aux_cut_graph, xk)
        elif model._settings["clique_method"] == "ostergard" :
            totalWeight, constr = clique.max_clique_ostergard(G, weights=xk, scale=10000)
        elif model._settings["clique_method"] == "ILS_VND" :
            totalWeight, constr = clique.max_clique_ILS_VND(G, weights=xk, scale=10000, iterations=10)
        if totalWeight > 1 + 1e-8 :
            model.cbLazy(grb.quicksum(model._x[u,v,k] for (u,v) in constr) <= 1)
    return


def bottom_up_no_cut_callback(model, where):
    if where == grb.GRB.Callback.MIPSOL :
        model.terminate()
    return


def check_settings(G, settings):
    """ Check the provided settings dictionary for errors and fill in defaults.

    Args:
        G (networkx.Graph):
            Simple connected graph.
        settings (dict):
            Dictionary containing settings options.

    Returns:
        None.

    Raises:
        ValueError:
            If the input is malformed.
    """
    defaults = {
        "K": G.order() - 1,
        "init_paths": dict(),  # i.e. all of them
        "add_symmtery_breaking": True,
        "GRB_symmtery_breaking": False,
        "GRB_verbose": True,
        "init_fix": "clique",
        "add_clique_cuts": True,
        "skip_dominated_pairs": True,
        "clique_method": "ostergard",
        "timeout": None,
        "print_solution": False,
        "num_heur_trials": None,
        "verbose": True,
    }
    for key in defaults:
        if key not in settings:
            settings[key] = defaults[key]
        elif key in settings and settings[key] is None:
            settings[key] = defaults[key]

    # Check settings here
    if settings["init_fix"].lower() not in ["lusp", "clique", "none"]:
        raise ValueError("'init_fix' must be 'lusp', 'clique', or 'none'")

    if settings["clique_method"].lower() not in ["ostergard", "brute", "ils_vnd", "none"]:
        raise ValueError("'clique_method' must be 'ostergard', 'brute', 'ILS_VND', or 'none'")


def graph_preprocessing(G, settings):
    """ Perform various preprocessing steps required for all solution methods.
    """
    # Build auxiliary cut graph, get skip pairs (among other things)
    prepr.presolve(G)

    # Build C/C++ graphs as necessary
    if settings["clique_method"] == "ILS_VND":
        G.Hcpp = clique.build_cppgraph(G)
        G.Hcpp.sort() # Is this necessary?
    elif settings["clique_method"] == "ostergard":
        G.Hc = clique.build_c_graph(G)

    # If you want to solve the path with an initial fixing
    if settings["init_paths"] == "heur":
        settings["init_paths"] = dict()
        for (u, v) in G.vertex_pairs:
            path = nx.shortest_path(G, source=u, target=v) # As nodes
            settings["init_paths"][u, v] = (path,)

    # Compute clique/lusp bound
    if settings["init_fix"] == "lusp":
        G.init_fix = prepr.find_longest_unique_shortest_path(G)
    elif settings["init_fix"] == "clique":
        if settings["clique_method"] == "brute":
            _, G.init_fix = clique.enumerate_best_clique(G.aux_cut_graph)
        elif settings["clique_method"] == "ostergard" :
            _, G.init_fix = clique.max_clique_ostergard(G, weights=None, scale=10000)
        elif settings["clique_method"] == "ILS_VND" :
            _, G.init_fix = clique.max_clique_ILS_VND(G, weights=None, scale=10000, iterations=10)


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


def model_preprocessing(G, model, settings):
    """ Fix the initial path """
    if settings["init_fix"].lower() != "none":
        fix_init_vars(G, model)
        if settings["verbose"] :
            print(f"Number of edges fixed by init_fix: {len(G.init_fix)}")


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
        help="Number of colors to use in the palatte (ignored if using bottom up)."
    )
    parser.add_argument(
        "-b", "--no-symmetry", action="store_true",
        help="Do not add color symmetry breaking constraints."
    )
    parser.add_argument(
        "-x", "--init-fix", type=str, default="clique",
        choices=["clique", "lusp"],
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

    settings = {
        "K": args.colors,
        "heur_bound": args.colors,
        "init_paths": None,  # None means add all shortest paths
        "add_symmetry_breaking": not args.no_symmetry,
        "GRB_symmetry_breaking": False,
        "init_fix": args.init_fix,
        "clique_method": args.clique_method,
        "skip_dominated_pairs": not args.no_skip_pairs,
        "add_clique_cuts": args.clique_cuts,
        "GRB_verbose": args.grb_verbose,
        "timeout": args.timeout,
        "print_solution": args.print_solution,
        "num_heur_trials": args.heur_trials,
        "verbose": args.verbose,
    }

    G = prepr.build_graph_from_file(args.filename)

    t0 = time.time()
    if args.method == "bu":
        srcg = src_bottom_up(G, settings)
    else:
        srcg = src_BnC(G, settings)
    if srcg is not None:
        print(f"src = {srcg:.1f}")
    else:
        print("Could not compute src")
    tf = time.time()
    print(f"Total time: {tf-t0:.2f}s")
