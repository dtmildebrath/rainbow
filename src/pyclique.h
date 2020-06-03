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
#ifndef PYCLIQUE_H
#define PYCLIQUE_H

#include "graph.h"

int maxclique_wrapper(int n, int m, int *edge_list, int *result);
graph_t* make_graph(int n, int m, int *edge_list);
void set_graph_weights(graph_t *g, int *weights);

#endif
