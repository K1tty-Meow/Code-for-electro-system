#!/bin/bash
#SBATCH -J MCMC
#SBATCH --nodelist=node7
#SBATCH -o mcmc_%A_%a.out    # File to which STDOUT will be written
#SBATCH -e mcmc_%A_%a.err    # File to which STDERR will be written
#SBATCH --cpus-per-task=1
#SBATCH -a 0-20

#c=$(echo "500.0"|bc)
input=($(cat sequ1.txt))
#input=($(seq 1.0 0.2 8.0))
ind=${input[$SLURM_ARRAY_TASK_ID]}
#g++ -O3 pmf.cpp -lm -o pmf.exe
./pmf.exe $ind
