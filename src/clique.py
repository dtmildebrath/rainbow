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
import ctypes

try:
    ccliquer = ctypes.CDLL("./libpyclique.so")
except OSError:
    print("Warning: no C cliquer library available")


## Brute force methods
def enumerate_best_clique(G, weights=None):
    """ Uses the networks find_cliques function to enumerate all cliques
        and simply take the largest clique. Returns the weight of the max weight clique
        and the vertices in that clique
        
        Should really only be for debugging
    """
    size = 0
    if weights is None:
        cliques = nx.algorithms.clique.find_cliques(G)
        best = None
        for clique in cliques:
            if len(clique) > size:
                size = len(clique)
                best = list(clique)
    else:
        cliques = nx.algorithms.clique.find_cliques(G)
        best = []
        for clique in cliques:
            if sum(weights[node] for node in clique) > size:
                best = list(clique)
                size = sum(weights[node] for node in clique)

    return size, best


def clique_stats(G):
    size, _ = enumerate_best_clique(G)
    counts = {i + 1: 0 for i in range(size)}
    cliques = nx.algorithms.clique.find_cliques(G)
    for clique in cliques:
        counts[len(clique)] += 1

    print(f"Graph of order {str(G.order())} has maximal cliques of sizes: ", counts)
    print(
        f"The aux cut graph has edge density {str(2*G.size()/(G.order()*(G.order()-1)))}"
    )
    return


## C version of Cliquer methods
def build_c_graph(G):
    set_bindings()
    H = G.aux_cut_graph

    m, n = H.size(), H.order()
    int_to_tuple = list(H.nodes())
    tuple_to_int = {val: i for (i, val) in enumerate(int_to_tuple)}
    G.int_to_tuple = int_to_tuple
    G.tuple_to_int = tuple_to_int

    edges = [[tuple_to_int[u], tuple_to_int[v]] for (u, v) in H.edges()]
    flattened_edges = [u for edge in edges for u in edge]

    c_n, c_m = ctypes.c_int(n), ctypes.c_int(m)
    c_edges = (ctypes.c_int * len(flattened_edges))(*flattened_edges)

    g = ccliquer.make_graph(c_n, c_m, c_edges)
    # g must be freed!!!!
    return g


class CGraph(ctypes.Structure):
    # Not clear how necessary this is
    _fields_ = [
        ("n", ctypes.c_int),
        ("weights", ctypes.POINTER(ctypes.c_int)),
    ]


def set_bindings():
    """ ccliquer is the return from ctypes.CDLL(path_name) """

    ## make_graph bindings
    ccliquer.make_graph.restype = ctypes.POINTER(CGraph)
    ccliquer.make_graph.argtypes = (
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(ctypes.c_int),
    )

    ## graph_free bindings
    ccliquer.graph_free.restype = None
    ccliquer.graph_free.argtypes = (ctypes.POINTER(CGraph),)

    ## maxclique_wrapper bindings
    ccliquer.maxclique_wrapper.restype = ctypes.c_int
    ccliquer.maxclique_wrapper.argtypes = (
        ctypes.POINTER(CGraph),
        ctypes.POINTER(ctypes.c_int),
    )


def max_clique_ostergard(G, weights=None, scale=10000):
    # Take care of the weights
    if weights is not None:
        weights = [int(scale * weights[tup]) for (i, tup) in enumerate(G.int_to_tuple)]
    else:  # Slightly slower but probably good to be careful
        weights = [1 for _ in range(len(G.int_to_tuple))]
        scale = 1

    c_weights = (ctypes.c_int * len(weights))(*weights)
    ccliquer.set_graph_weights(G.Hc, c_weights)

    # Set up the return values
    result = [0 for _ in range(G.Hc.contents.n + 1)]
    c_result = (ctypes.c_int * len(result))(*result)

    # Call the actual function
    weight = ccliquer.maxclique_wrapper(G.Hc, c_result)

    size = c_result[0]
    nodes = [G.int_to_tuple[v] for v in c_result[1 : size + 1]]
    return weight / scale, nodes


def free_graph(g):
    ccliquer.graph_free(g)


if __name__ == "__main__":
    import argparse
    import prepr
    import time
    parser = argparse.ArgumentParser(
        description="Compute the clique lower bound for src(G) (using Cliquer)"
    )

    parser.add_argument("filename", help="Path to DIMACS file")
    args = parser.parse_args()

    G = prepr.build_graph_from_file(args.filename)
    t0 = time.time()
    prepr.presolve(G) 
    t1 = time.time()
    print(f"Presolve time: {t1-t0:.4f}s")
    G.Hc = build_c_graph(G)
    bound, _ = max_clique_ostergard(G)
    t2 = time.time()
    print(f"Cliquer time: {t2-t1:.4f}s")
    print(f"Clique bound: {bound:.1f}")
    free_graph(G.Hc)
