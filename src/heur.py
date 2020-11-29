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
import random
import prepr

def path_fixing_heuristic(G, num_trials=3, online_sampling=False, fmethod="clique"):
    """ Compute a valid strong rainbow coloring of G.

    Args:
        G (Graph):
            Graph to compute a coloring for. Must have the following attached
            attributes:
            - `vertex_pairs`
            - `ordered_edges`
            - `diameter`
        num_trials (int, optional):
            Number of trials to execute. Default is 3.
        online_sampling (boolean, optional):
            If True, randomly samples a shortest path between each pair of
            vertices as needed. Otherwise, enumerate every shortest path in G
            before starting the heuristic. Default is False.
        fmethod (str, optional):
            Method used to fix an initial set of edges. Must be either `clique`
            or `path`. If `clique`, a set of edges corresponding to a clique in
            the auxiliary graph H are fixed. Otherwise, a set of edges in
            longest shortest path in G are fixed (i.e. a shortest path whose
            length equals diam(G)). Default is `clique`.

    Returns:
        dict:
            Dictionary mapping edges to colors (integers), corresponding to a
            valid strong rainbow coloring of G.
        int:
            The number of colors used in the strong rainbow coloring.
    """
    required_attrs = ["vertex_pairs", "ordered_edges", "diameter"]
    for attr in required_attrs:
        if not hasattr(G, attr):
            raise RuntimeError(
                f"Graph must have a '{attr}' attribute before running path fixing heuristic"
            )
    if fmethod not in ("clique", "path"):
        raise ValueError("fmethod must `clique` or `path`")

    best_bound = G.size()+1
    best_paths = {pair:0 for pair in G.vertex_pairs}
    best_colors = {e:0 for e in G.ordered_edges}
    diam = G.diameter

    if not online_sampling:
        # Enumerate all shortest paths using
        # brute force so that we can sample them
        # This is quite a bit faster than explicitly
        # computing and storing all of the pairwise DAGS
        # (and sampling is simpler + faster)
        # This takes <1s for ieee118
        shortest_paths = dict()
        for (u, v) in G.vertex_pairs:
            path_list = list()
            paths = nx.all_shortest_paths(G, source=u, target=v)
            for path in paths:
                path_list.append(
                    prepr.path_to_edge_tuples(path, G.ordered_edges) 
                )
            shortest_paths[u, v] = path_list

    # Run the heuristic
    for trial in range(num_trials):
        color_fixings = {e: -1 for e in G.ordered_edges}
        path_fixing = {pair: None for pair in G.vertex_pairs}
        K = 0

        # Choose a path fixing
        if online_sampling:
            for (u,v) in G.vertex_pairs:
                path_fixing[u, v] = _sample_shortest_path(G, u, v)
        else:
            for (u,v) in G.vertex_pairs:
                path_ix = random.randint(0, len(shortest_paths[u,v])-1)
                path_fixing[u,v] = shortest_paths[u,v][path_ix]

        if fmethod == "path":
            # Fix colors along some longest path
            for (u,v) in G.vertex_pairs:
                if len(path_fixing[u,v]) == diam:
                    for edge in path_fixing[u,v]:
                        color_fixings[edge] = K
                        K += 1
                    break
        elif fmethod == "clique":
            # Fix the edges in the maximum clique in H
            for edge in G.init_fix:
                color_fixings[edge] = K
                K += 1

        free_edges = [edge for edge in G.ordered_edges if color_fixings[edge] == -1]
        random.shuffle(free_edges)
        
        # K is the number of colors we have used so far
        for free_edge in free_edges:
            color_available = [True for _ in range(K)]
            num_available_colors = K
            for (u, v) in G.vertex_pairs:
                if num_available_colors == 0:
                    break
                path = path_fixing[u,v]
                if free_edge in path:
                    for edge in path:
                        color = color_fixings[edge]
                        if color != -1 and color_available[color]:
                            color_available[color] = False
                            num_available_colors -= 1
            if num_available_colors > 0:
                options = [k for (k, av) in enumerate(color_available) if av]
                color_ix = random.randint(0, len(options)-1)
                new_color = options[color_ix]
            else:
                new_color = K
                K += 1
            color_fixings[free_edge] = new_color

        # After all edges are colored, save solution if it beats best_bound
        if K < best_bound:
            best_bound = K
            best_sol = color_fixings
            best_paths = path_fixing

            if K == diam:
                print("Heuristic found an optimal solution!")
                break
                # return best_sol, best_bound

    # print(best_sol)
    return best_sol, best_bound


def _sample_shortest_path(G, source, target):
    """ Sample a path uniformly at random from the set of shortest paths
    connecting source and target in G.

    Args:
        G (Graph):
            Networkx graph.
        source (object):
            Source node. Must be an element of the node set of G.
        target (object):
            Target node. Must be an element of the node set of G.

    Returns:
        list of tuple:
            Resulting path stored as a list of edge tuples.

    Notes:
    - Each node in the graph G must have an attached directed acyclic graph
      (DAG), stored with the key `dag`. The DAG must have an integer `r` stored
      at each node, corresponding to the number of shortest paths in G.
    """
    D = G.node[source]["dag"]
    path = list()
    L = tuple(D.predecessors(target))
    weights = [D.node[u]["r"] for u in L]
    node = random.choices(L, weights=weights)[0]
    path = [node, target]
    while node != source:
        L = tuple(D.predecessors(node))
        weights = [D.node[u]["r"] for u in L]
        node = random.choices(L, weights=weights)[0]
        path.insert(0, node)
    path = prepr.path_to_edge_tuples(path, G.ordered_edges) 
    return path


if __name__ == "__main__":
    import time
    import src
    import argparse
    import math

    parser = argparse.ArgumentParser(
        description="Produce a valid strong rainbow coloring of G, which is an upper bound on src(G)"
    )

    parser.add_argument("filename", help="Path to DIMACS file")
    parser.add_argument(
        "-p", "--print-solution", action="store_true",
        help="Display the resulting edge coloring"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Turn on output",
    )
    parser.add_argument(
        "-x", "--init-fix", type=str, default="clique",
        choices=["clique", "path"],
        help="Type of initial edge fixing to use. Default is `clique`."
    )
    parser.add_argument(
        "-n", "--trials", type=int,
        help="Number of heuristic trials to run. "
        "Default is the number of vertices in the graph divided by 5."
    )
    parser.add_argument(
        "-o", "--online-sampling", action="store_true",
        help="Use online sampling (offline is faster as number of trials increases)"
    )

    args = parser.parse_args()

    settings = {
        "clique_method": "ostergard" if args.init_fix == "clique" else None,
        "init_fix": args.init_fix,
        "init_paths": None,
    }

    G = prepr.build_graph_from_file(args.filename)
    if args.verbose:
        print("Starting preprocessing...", end="", flush=True)
    src.graph_preprocessing(G, settings)
    if args.verbose:
        print("done.")

    if args.trials is None:
        args.trials = math.ceil(G.order() / 5)

    start = time.time()
    cmap, rc = path_fixing_heuristic(
        G, 
        num_trials=args.trials,
        fmethod=args.init_fix,
        online_sampling=args.online_sampling,
    )
    end = time.time()
    print(f"Bound: {rc}")
    print(f"Time: {end-start:.4f}s")
    if args.print_solution:
        for (u, v) in cmap:
            print(f"Edge ({u:d},{v:d}) is color {cmap[u,v]:d}")
