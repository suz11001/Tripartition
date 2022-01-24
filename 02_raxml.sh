#!/bin/bash
#SBATCH --job-name=rax
#SBATCH --mail-user=sumaira.zaman@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=8G
#SBATCH -o raxml_%j.out
#SBATCH -e raxml_%j.err
#SBATCH --partition=general
#SBATCH --qos=general

module load anaconda
module load RAxML/8.2.11

for w in $1/window*; do
    ls -lh $w
    window=$(basename $w)
    echo $window
    raxmlHPC-PTHREADS -T 6 -f a -x 12345 -p 12345 -s $w -n $window -m $2 -#100 
done

