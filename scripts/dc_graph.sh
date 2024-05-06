#!/bin/bash --login

#SBATCH --time=04:00:00   # walltime
#SBATCH --output=comp_results_sage.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=4096M   # memory per CPU core
#SBATCH -J "simple_comp_graph_test"   # job name
#SBATCH --mail-user=dseyfr99@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
cd /home/seyfdall/compute/graph_theory/laplace-spread-conjecture
mamba activate sage_mpi
# module load mpi/openmpi-1.10.7_gnu4.8.5
# export MPICC=$(which mpicc)
python3 dc_graph_comparison.py