/*  This file is part of SRC.
 *
 *  SRC is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SRC is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SRC.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Copyright 2020 David T. Mildebrath, Logan A. Smith
 */
#include <stdio.h>
#include "cliquer.h"

void set_graph_weights(graph_t *g, int *weights)
{
    int i;
    for (i=0; i<g->n; i++) {
        g->weights[i] = weights[i];
    }
}

// There is already a graph_free method

graph_t* make_graph(int n, int m, int *edge_list)
{
    graph_t *g;
    int i;

    g = calloc(1, sizeof(graph_t));
    g->n = n;

    // Set up the empty graph with the right number of nodes
    g->edges = calloc(g->n, sizeof(set_t)); // Each node has list of adjacent nodes
    for (i=0; i<g->n; i++) {
        g->edges[i] = set_new(g->n); // Each node can have at most n neighbors
    }
    g->weights = calloc(g->n, sizeof(int)); // Integer weights
    for (i=0; i<g->n; i++) {
        g->weights[i] = 1;
    }

    // Add the edges
    for (i=0; i<2*m; i+=2) {
        GRAPH_ADD_EDGE(g, edge_list[i], edge_list[i+1]); 
    }

    // Check the graph
    ASSERT(graph_test(g, NULL)); // Could point to stderr not NULL

    return g;
}

int maxclique_wrapper(graph_t *g, int *result)
{
    if (result == NULL)
        fprintf(stderr, "WE GOOFED LEAKS EVERYWHERE\n");
    int i, j, omega, weight; // omega = size (not weight!) of max clique
    set_t s;

    // Find the max clique
    clique_default_options->time_function = NULL;
    s = clique_find_single(g, 0, 0, FALSE, clique_default_options);

    omega = set_size(s);
    weight = 0;
    result[0] = omega;
    j = 1;
    for (i=0; i<SET_MAX_SIZE(s); i++) {
        if (SET_CONTAINS(s, i)) {
            result[j] = i;
            weight += g->weights[i];
            j++;
        }
    }

    set_free(s);
    return weight;
}
