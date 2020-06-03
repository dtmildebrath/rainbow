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
import src
import prepr
import time


def run_tests() :
    path = "../graphs"
    heur_bounds = _load_heur_file("heur_bounds.csv")
    inst_list = prepr.read_instance_names("insts.txt")

    test_bottom_up_cuts    = 0
    test_bottom_up_no_cuts = 0
    test_BnC_cuts          = 0
    test_BnC_no_cuts       = 0
    test_BnC_naive         = 1

    time_limit = 3600
    graph_statistics = 1
    verbose = 1

    times = {inst:{} for inst in inst_list}

    for inst in inst_list :
        print(f"\n##############################################################\n")
        print(f"Beginning tests for instance: {inst}")
        G = prepr.build_graph_from_file(f"{path}/{inst}")
        src_val = None
        heur_bound = heur_bounds[inst]["bound"]
        heur_time = heur_bounds[inst]["time"]

        if test_bottom_up_cuts :
            settings = {
                "K": None,
                "init_paths": None,  # None means add all shortest paths
                "add_symmetry_breaking": True,
                "GRB_symmetry_breaking": False,
                "init_fix": "clique",
                "clique_method": "ostergard",
                "skip_dominated_pairs": True,
                "add_clique_cuts": True,
                "GRB_verbose": verbose,
                "timeout": time_limit
            }
            start = time.time()
            sol = src.src_bottom_up(G, settings)
            end = time.time()
            times[inst]["bottom_up_cuts"] = end-start
            if sol is not None :
                if src_val is not None and sol != src_val :
                    print(f"METHODS RETURNED DIFFERENT SOLUTIONS: {src_val} and {sol}")
                src_val = sol
            if sol is None :
                print("--- Method: bottom_up_cuts Timed Out!")

        if test_bottom_up_no_cuts :
            settings = {
                "K": None,
                "init_paths": None,  # None means add all shortest paths
                "add_symmetry_breaking": True,
                "GRB_symmetry_breaking": False,
                "init_fix": "clique",
                "clique_method": "ostergard",
                "skip_dominated_pairs": True,
                "add_clique_cuts": False,
                "GRB_verbose": verbose,
                "timeout": time_limit
            }
            start = time.time()
            sol = src.src_bottom_up(G, settings)
            end = time.time()
            times[inst]["bottom_up_no_cuts"] = end-start
            if sol is not None :
                if src_val is not None and sol != src_val :
                    print(f"METHODS RETURNED DIFFERENT SOLUTIONS: {src_val} and {sol}")
                src_val = sol
            if sol is None :
                print("--- Method: bottom_up_no_cuts Timed Out!")

        if test_BnC_cuts :
            settings = {
                "K": None,
                "init_paths": None,  # None means add all shortest paths
                "add_symmetry_breaking": True,
                "GRB_symmetry_breaking": False,
                "init_fix": "clique",
                "clique_method": "ostergard",
                "skip_dominated_pairs": True,
                "add_clique_cuts": True,
                "GRB_verbose": verbose,
                "timeout": time_limit,
                "heur_bound": heur_bound,
            }
            start = time.time()
            sol = src.src_BnC(G, settings)
            end = time.time()
            times[inst]["BnC_cuts"] = end-start+heur_time
            if sol is not None :
                if src_val is not None and sol != src_val :
                    print(f"METHODS RETURNED DIFFERENT SOLUTIONS: {src_val} and {sol}")
                src_val = sol
            if sol is None :
                print("--- Method: BnC_cuts Timed Out!")


        if test_BnC_no_cuts :
            settings = {
                "K": None,
                "init_paths": None,  # None means add all shortest paths
                "add_symmetry_breaking": True,
                "GRB_symmetry_breaking": False,
                "init_fix": "clique",
                "clique_method": "ostergard",
                "skip_dominated_pairs": True,
                "add_clique_cuts": False,
                "GRB_verbose": verbose,
                "timeout": time_limit,
                "heur_bound": heur_bound,
            }
            start = time.time()
            sol = src.src_BnC(G, settings)
            end = time.time()
            times[inst]["BnC_no_cuts"] = end-start+heur_time
            if sol is not None :
                if src_val is not None and sol != src_val :
                    print(f"METHODS RETURNED DIFFERENT SOLUTIONS: {src_val} and {sol}")
                src_val = sol
            if sol is None :
                print("--- Method: BnC_no_cuts Timed Out!")

        if test_BnC_naive :
            settings = {
                "K": None,
                "init_paths": None,  # None means add all shortest paths
                "add_symmetry_breaking": True,
                "GRB_symmetry_breaking": False,
                "init_fix": "none",
                "clique_method": "none",
                "skip_dominated_pairs": False,
                "add_clique_cuts": False,
                "GRB_verbose": verbose,
                "timeout": time_limit,
                "heur_bound": heur_bound,
            }
            start = time.time()
            sol = src.src_BnC(G, settings)
            end = time.time()
            times[inst]["BnC_naive"] = end-start+heur_time
            if sol is not None :
                if src_val is not None and sol != src_val :
                    print(f"METHODS RETURNED DIFFERENT SOLUTIONS: {src_val} and {sol}")
                src_val = sol
            if sol is None :
                print("--- Method: BnC_naive Timed Out!")

        if graph_statistics :
            src.graph_statistics(G)
            if src_val is not None :
                print(f"Graph src: {int(src_val)}")

        print(f"\nTiming Results for {inst}:")
        for method in times[inst] :
            print(f"\t {method}: \t {times[inst][method]:.4f}")
        print(f"\n##############################################################\n\n")


def print_results(test_insts, times) :
    for inst in test_insts :
        print(f"\nInstance {inst}:")
        for method in times[inst] :
            print(f"\t {method}: \t {times[inst][method]:.4f}")


def _load_heur_file(fname):
    bounds = dict()
    try:
        with open(fname, "r") as f:
            f.readline()
            for line in f:
                line = line.split(",")
                bounds[line[0]] = {
                    "bound": int(line[1]),
                    "trials": int(line[2]),
                    "time": float(line[3]),
                }   
    except FileNotFoundError:
        print(
            "Must create a heur_bounds.csv file in this directory.\n"
            "Do this by running this command:\n\n"
            "    python preheur.py >> heur_bounds.csv\n\n"
            "This command may take several minutes to execute."
        )
        quit()
    return bounds


if __name__ == "__main__" :
    run_tests()
