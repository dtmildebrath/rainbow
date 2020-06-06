Copyright 2020 David T. Mildebrath, Logan A. Smith

This code is accomponied by [this paper](https://arxiv.org/abs/2006.02988).

# Strong Rainbow Connection (SRC)

This code provides several methods for computing the strong rainbow connection
numbers of general (simple) graphs. This code is available under the GNU
General Public License v3.0.

## Installation Instructions

A large majority of this code relies on the
[Gurobi solver](https://www.gurobi.com "Gurobi solver"). Free versions of this
solver are available for academics.

The code also relies on the open source
[Cliquer](https://users.aalto.fi/~pat/cliquer.html "Cliquer") code. As our code
is written in the Python programming language, and Cliquer is written in the C
programming language, we provide some simple Python bindings for Cliquer. Steps
to use our bindings are as follows:

1. Download and unzip the
   [Cliquer tar.gz file](http://users.aalto.fi/~pat/cliquer/cliquer-1.21.tar.gz "Download Cliquer").
2. Move the files `pyclique.c` and `pyclique.h` that we provide into the
   directory with the Cliquer source.
3. Add the following lines to the Cliquer Makefile:
   
        pyclique: pyclique.o cliquer.o graph.o reorder.o
            gcc -fPIC -c -Wall pyclique.c cliquer.c graph.c reorder.c
            ld -shared pyclique.o cliquer.o graph.o reorder.o -o libpyclique.so

4. Run `make pyclique`
5. Move the resulting `libpyclique.so` file to the directory containing the
   Python source code (i.e. `prepr.py`, `src.py`, etc.).


## Usage Instructions

We provide three separate command line utilities for computing src(G). All
methods assume that the specified graph is given in DIMACS format.

1. `src.py` exactly computes the src of a graph, using either branch-and-cut or
   bottom-up methods. Type `python src.py -h` for full details.
2. `heur.py` computes a heuristic solution (i.e. upper bound) for the src(G).
   Type `python heur.py -h` for full details.
3. `clique.py` computes the clique lower bound w'(G) for src(G). This method
   uses Cliquer for exactly computing the maximum clique in the auxiliary graph
   H. Type `python clique.py -h` for full details.

## Reproducing Test Results from the Paper

The script `test.py` reproduces the test results given in the paper. In order
to ensure fair testing, we only run the randomized heuristic once for each
instance and save the results, then use the saved results for all tests.
Consequently, before running the `test.py` script, run the command

        python preheur.py >> heur_bounds.csv

to execute the heuristic and save the results. To change which instances are
tested, modify `insts.txt` appropriately. Which tests are run can be modified
by changing the 0/1 flags in the `run_tests()` function of `test.py`
appropriately.
