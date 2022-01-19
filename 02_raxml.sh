#!/bin/bash
#SBATCH --job-name=rax-de-w
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
#SBATCH --array=1-100%35

module load anaconda
module load RAxML/8.2.11

i=$SLURM_ARRAY_TASK_ID


d="/home/FCAM/szaman/researchMukul/pgtr/nov8_2020/window_analysis_60perc/split/sequences_control_domain60perc/"
echo $d
dir=$(basename $d)
echo $dir
for paramPath in $d/*/; do 
    echo $paramPath
    param=$(basename $paramPath)
    mkdir -p ~/researchMukul/pgtr/nov8_2020/window_analysis_60perc/trees/$dir/$param/$i
    cd ~/researchMukul/pgtr/nov8_2020/window_analysis_60perc/trees/$dir/$param/$i
    for w in ~/researchMukul/pgtr/nov8_2020/window_analysis_60perc/split/$dir/$param/$i/window*; do
	ls -lh $w
	window=$(basename $w)
	echo $window
	mkdir -p $window
	cd $window
	echo "I am here:"
	echo $PWD
	raxmlHPC-PTHREADS -T 6 -f a -p 13579 -N 100 -m GTRCAT -x 12345 -s $w -n $window
	cd ~/researchMukul/pgtr/nov8_2020/window_analysis_60perc/trees/$dir/$param/$i
    done
    cd ~/researchMukul/pgtr/nov8_2020/window_analysis_60perc/trees/$dir/
done

