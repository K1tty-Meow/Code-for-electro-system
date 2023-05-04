#!/bin/bash
#SBATCH -J MCMC
#SBATCH --nodelist=node5
#SBATCH -o mcmc_%A_%a.out    # File to which STDOUT will be written
#SBATCH -e mcmc_%A_%a.err    # File to which STDERR will be written
#SBATCH --cpus-per-task=1
#SBATCH -a 0-21

#c=$(echo "500.0"|bc)
c=100000
input=($(cat sequ2.txt))
#input=($(seq 1.0 0.2 8.0))
ind=${input[$SLURM_ARRAY_TASK_ID]}
#g++ -O3 numbernew.cpp -lm -o newnu.exe
./newnu.exe $ind $c
result=`python sample.py $ind`
while (($result > 0))
do
#	g++ -O3 numbernew.cpp -lm -o newnu.exe
	let c+=5000
	./newnu.exe $ind $c
	result=`python sample.py $ind`
done

