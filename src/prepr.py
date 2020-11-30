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
""" General utility functions for the SRC problem
"""
import networkx as nx
import itertools

from collections import deque

def presolve(G, skip_dominated_pairs=True):
    """ Compute DAGs, auxiliary graph, and other graph parameters used when
    computing src(G).

    This function does a lot of things, because many of them are tightly
    coupled.

    Args:
        G (Graph):
            A simple, connected graph.
        skip_dominated_pairs (bool, optional):
            If True, G.vertex_pairs does not include pairs of vertices which do
            not need to be explicitly enforced to be rainbow connected in the
            IP model. Default is True.
    
    Returns:
        None

    Attaches to G:
        ordered_edges (list of tuple):
            List of all edges (u, v) in E(G) with u < v.
        path_counts (dict):
            A map from non-adjacent vertex pairs (u, v) to the number of
            shortest paths between u and v in G.
        aux_cut_graph (Graph):
            The auxiliary graph H (whose nodes correspond to edges of G).
        vertex_pairs (tuple of tuple):
            The pairs of vertices which must be explicitly considered when
            computing src(G) (i.e. vertex pairs after removals).
        diameter (int):
            The diameter of the graph G.

    Attaches to each node u in G:
        dag (DiGraph):
            The directed graph rooted at node u.
        r (dict):
            G.node[u]["r"][v] is the number of shortest u-v paths in G. 
    """
    # Build the (single source) DAGs
    _dag_prep(G)
    G.ordered_edges = get_ordered_edges(G)

    # All non-adjacent vertex pairs (u, v) in G
    # Sorted so that u < v
    v_pairs = tuple(
        (min(u, v), max(u, v))
        for (u, v) in itertools.combinations(G.nodes(), 2)
        if u not in G.neighbors(v)
    )

    # Initialize data structures
    H_edges = dict() # Really just a hash table, values don't matter
    skip_pairs = dict() # Also just a hash table
    path_counts = {(u, v): 0 for (u, v) in v_pairs}

    # Iterate over pairs and compute edges
    for (source, sink) in v_pairs:
        path_count, cut_edges, cut_vertices = _presolve_helper(G, source, sink)
        path_counts[source, sink] = path_count

        if len(cut_edges) >= 2:
            for (u, v) in itertools.combinations(cut_edges, 2):
                H_edges[u,v] = True # Value does not matter

        if skip_dominated_pairs:
            for v in cut_vertices:
                if source not in G.neighbors(v):
                    skip_pairs[min(source, v), max(source, v)] = True
                if sink not in G.neighbors(v):
                    skip_pairs[min(sink, v), max(sink, v)] = True

    remaining_pairs = tuple(pair for pair in v_pairs if pair not in skip_pairs)

    # Build graph object
    H = nx.Graph()
    H.add_nodes_from(G.ordered_edges)
    H.add_edges_from(H_edges.keys())

    # Attach the results (do not return)
    G.path_counts = path_counts
    G.aux_cut_graph = H
    G.vertex_pairs = remaining_pairs
    G.diameter = nx.diameter(G)


def _dag_prep(G):
    """ Construct and store all the DAGs
    """
    for v in G.nodes():
        G.node[v]["r"] = {node: 0 for node in G.nodes()}

    for v in G.nodes():
        D = nx.DiGraph()
        D.add_edges_from(_build_dag_edges(G, v))
        G.node[v]["dag"] = D

    # Copy the path counts into the Dags
    # If D = G.node[v]["d"]
    # Then D.node[u]["r"] is the number of shortest u, v paths in G
    for source in G.nodes():
        D = G.node[source]["dag"]
        for u in G.nodes():
            D.node[u]["r"] = G.node[source]["r"][u]


def _build_dag_edges(G, source, neighbors=None, depth_limit=None):
    """ Takes a graph and a vertex source. Returns an iterator over the edges
        in the directed graph construction rooted at source. 

        Each edge has a weight, corresponding to the distance of
        that edge from the source node (edges incident upon the source node
        have distance 1, etc.). These numbers are used later when computing H.

        Commented out yield statements return 1-based layer num. of current edge

        Populates a dictionary at stored at the source node in G (not D!)
        G.node[source]["r"] is a dict which maps nodes to ints.
        D.node[source]["r"][v] equals the number of distinct shortest source, v
        paths in G.
    """
    visited = {source}
    G.node[source]["r"][source] = 1
    if depth_limit is None:
        depth_limit = len(G)
    depthOf = {source: depth_limit}
    queue = deque([(source, depth_limit, G.neighbors(source))])
    while queue:
        parent, depth_now, children = queue[0]
        try:
            child = next(children)
            if child in visited and depthOf[child] == depth_now - 1:
                G.node[source]["r"][child] += G.node[source]["r"][parent]
                yield parent, child # , depth_limit - depth_now + 1

            if child not in visited:
                G.node[source]["r"][child] = G.node[source]["r"][parent]
                yield parent, child # , depth_limit - depth_now + 1
                visited.add(child)
                depthOf[child] = depth_now - 1
                if depth_now > 1:
                    queue.append((child, depth_now - 1, G.neighbors(child)))
        except StopIteration:
            queue.popleft()


def _presolve_helper(G, source, sink):
    """ This function does three things (as all good functions do):
        
        1. Count the number of distinct shortest paths between source and sink
        2. Compute the "cut edges" between source and sink, i.e., the edges in
            G used by all shortest paths between source and sink
        3. Compute the "cut vertices" between source and sink, i.e., the
            vertices in G used by all shortest paths between source and sink
    """
    D = G.node[source]["dag"]
    L = set([sink])
    D.node[sink]["r"] = 1
    cut_edges = list()
    cut_vertices = list()

    while source not in L:
        N = set(x for j in L for x in D.predecessors(j))
        for j in N:
            D.node[j]["r"] = sum(D.node[v]["r"] for v in D.successors(j) if v in L)
        if len(N) == 1 and len(L) == 1:
            u, v = next(iter(L)), next(iter(N))
            cut_edges.append((min(u, v), max(u, v)))
        if len(N) == 1 and source not in N:
            cut_vertices.append(next(iter(N)))
        L = N
    return D.node[source]["r"], cut_edges, cut_vertices


def path_to_edge_tuples(path, full_edge_list):
    """ Compute a list of edge tuples contained in a 
        path given a list of nodes contained in that path. 

        Has an extra check because the ordering of the edge tuples 
        in the original (undirected) graphs is essentially random.
    """
    edges = [(path[i], path[i + 1]) for i in range(len(path) - 1)]
    for (ix, (u, v)) in enumerate(edges):
        if (u, v) not in full_edge_list:
            edges[ix] = (v, u)
    return edges


def count_all_paths(G):
    """ Compute the total number of shortest paths between all unique vertex pairs
        in G, and return the result as a dictionary.

        This function is not necessary for pre-processing; it is
        only used for data collection. It is not optimized for speed.

        This function does not require any preprocessing of the graph G
        whatsoever (not even G.vertex_pairs).
    """
    nodes = list(G.nodes())
    path_counts = dict()
    for i in range(len(nodes)):
        for j in range(i+1, len(nodes)):
            u, v = min(nodes[i], nodes[j]), max(nodes[i], nodes[j])
            path_counts[u, v] = len(list(
                nx.all_shortest_paths(G, source=u, target=v)
            ))
    return path_counts


def get_ordered_edges(G):
    return [(u, v) if u < v else (v, u) for (u, v) in G.edges()]


def find_longest_unique_shortest_path(G):
    """ Hell of a function name, huh.
    """
    longest_path = None
    diameter = 0
    for (u, v) in G.vertex_pairs:
        num_paths = G.path_counts[u, v]
        if num_paths == 1:
            path = nx.shortest_path(G, source=u, target=v)
            if len(path) > diameter:
                diameter = len(path)
                longest_path = path_to_edge_tuples(path, G.ordered_edges)
    return longest_path


def build_graph_from_file(fname):
    if fname[-4:] == '.gml' :
        G = nx.read_gml(fname, label='id')

    elif fname[-4:] == '.txt' :
        with open(fname, "r") as f:
            m = int(f.readline().strip().split()[-1])
            edges = [None for _ in range(m)]
            for index in range(m):
                line = f.readline().strip().split()
                u, v = int(line[0]), int(line[1])
                edges[index] = (u, v)
        G = nx.Graph()
        G.add_edges_from(edges)

    else:
        raise ValueError(f"File extension of {fname} not recognized, unsure of how to build graph.")

    return G


def read_instance_names(fname):
    """ Read a plain text file of graph instance names (including file
        extensions), one instance per line.
    """
    insts = list()
    with open(fname, "r") as f:
        for line in f:
            insts.append(line.strip())
    return insts

if __name__ == "__main__":
    print("prepr.py has no main")
