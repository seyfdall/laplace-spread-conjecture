#!/bin/bash --login

#SBATCH --time=00:10:00   # walltime
#SBATCH --output=mpi_comp_results_sage.txt
#SBATCH --nodes=1
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=2048M   # memory per CPU core
#SBATCH -J "mpitest"   # job name
#SBATCH --qos=test

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
cd /home/seyfdall/compute/graph_theory/laplace-spread-conjecture
mamba activate sage_mpi
# module load mpi/openmpi-1.10.7_gnu4.8.5
export MPICC=$(which mpicc)
mpirun -np 2 python3 dc_graph_comparison.py