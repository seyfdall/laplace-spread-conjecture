# Laplace-Spread-Conjecture
This Repository is designed to numerically test and visualize ideas for alternate paths to proving the Laplacian Spread Conjecture and stronger conjectures as well 
(see https://arxiv.org/abs/2201.04225 for the motivation behind this project).  It utilizes Sagemath based in Python for graph generation and MPI to parallelize numerical tests.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
- [Contributing](#contributing)
- [License](#license)
- [Credits](#credits)

## Installation
Steps for downloading the Github Repo and setting up a working mamba environment.  As of writing, MPI and Sage were not compatible over Python 3.12, so the more stable version
3.11 was used.

```bash
$ git clone https://https://github.com/seyfdall/laplace-spread-conjecture.git
$ cd laplace-spread-conjecture
$ mamba create -n sage_mpi python=3.11
$ mamba activate sage_mpi
$ mamba install openmpi
$ mamba install mpi4py
$ mamba install sage
```

## Usage
To test your setup, run the following command after activating your mamba environment.

```bash
$ sbatch mpitest.sh
```

This should output two lines of Hello World to a file along with what process number they are.  If this works, MPI is up and running and you can move onto the sage code.
Since the sage code is packaged into your mamba environment, you should also be able to run the comp_graph_comparison.py file which as of writing will print the dominating
edges and output a png image into your directory.

```bash
$ python comp_graph_comparison.py
```
