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
import heur
import prepr
import time
import src
import clique
import math


def main():
    path = "../graphs"
    inst_list = prepr.read_instance_names("insts.txt")

    settings = {
        "clique_method": "ostergard",
        "init_fix": "clique",
        "init_paths": None,
    }

    print("Instance,Bound,Trials,Time(s)")
    for fname in inst_list:
        G = prepr.build_graph_from_file(f"{path}/{fname}")
        src.graph_preprocessing(G, settings)
        num_trials = math.ceil(G.order() / 5)
        t0 = time.time()
        _, rc = heur.path_fixing_heuristic(G, num_trials=num_trials)
        t1 = time.time()
        print(f"{fname},{rc},{num_trials},{t1-t0}")
        clique.free_graph(G.Hc)


if __name__=="__main__":
    main()
