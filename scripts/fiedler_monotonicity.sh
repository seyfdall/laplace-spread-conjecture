#!/bin/bash --login

#SBATCH --time=00:20:00   # walltime
#SBATCH --output=./fiedler_monotonicity_results/results.txt
#SBATCH --nodes=1
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=1024M   # memory per CPU core
#SBATCH -J "fiedler_montonicity"   # job name
#SBATCH --mail-user=dseyfr99@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
cd /nobackup/autodelete/usr/seyfdall/graph_theory/laplace-spread-conjecture
mamba activate sage_mpi
# module load mpi/openmpi-1.10.7_gnu4.8.5
# export MPICC=$(which mpicc)
python3 fiedler_monotonicity.py