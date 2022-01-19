#!/bin/bash
#SBATCH --job-name=rax-bs
#SBATCH --mail-user=sumaira.zaman@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=8G
#SBATCH -o raxml_bs%j.out
#SBATCH -e raxml_bs%j.err
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=1-100%50

module load anaconda
module load RAxML/8.2.11

i=$SLURM_ARRAY_TASK_ID

# for d in ~/researchMukul/pgtr/nov8_2020/window_analysis/01_split/sequences_control_domain*perc/; do
#     echo $d
#     dir=$(basename $d)
#     echo $dir
#     if  [[ $dir == *"sequences_control_domain40perc"* ]]; then
# 	for paramPath in $d/*; do 
# 	    param=$(basename $paramPath)
# 	    if [[ $param == *"control"* ]]; then
# 		mkdir -p ~/researchMukul/pgtr/nov8_2020/window_analysis/03_bootstrap_trees/$dir/$param/$i/
# 		cd ~/researchMukul/pgtr/nov8_2020/window_analysis/03_bootstrap_trees/$dir/$param/$i/
# 		for w in ~/researchMukul/pgtr/nov8_2020/window_analysis/01_split/$dir/$param/$i/window*; do
# 		    if  [[ $w != *"reduced"* ]]; then
# 			ls -lh $w
# 			window=$(basename $w)
# 			echo $window
# 			mkdir -p $window
# 			cd $window
# 			echo "I am here:"
# 			echo $PWD
# 			raxmlHPC-PTHREADS -T 6 -b 12345 -p 13579 -# 100 -m GTRCAT -s $w -n $window
# 			cd ~/researchMukul/pgtr/nov8_2020/window_analysis/03_bootstrap_trees/$dir/$param/$i
# 		    fi
# 		done
# 		break 3
# 	    fi
# 	done
#     fi
# done


for d in ~/researchMukul/pgtr/nov8_2020/window_analysis_60perc/split/sequences_control_domain*perc/; do
    echo $d
    dir=$(basename $d)
    echo $dir
    if  [[ $dir == *"sequences_control_domain60perc"* ]]; then
        for paramPath in $d/*; do
            param=$(basename $paramPath)
            if [[ $param == *"control"* ]]; then
                mkdir -p ~/researchMukul/pgtr/nov8_2020/window_analysis_60perc/03_bootstrap_trees/$dir/$param/$i/
                cd ~/researchMukul/pgtr/nov8_2020/window_analysis_60perc/03_bootstrap_trees/$dir/$param/$i/
                for w in ~/researchMukul/pgtr/nov8_2020/window_analysis_60perc/split/$dir/$param/$i/window*; do
                    if  [[ $w != *"reduced"* ]]; then
                        ls -lh $w
                        window=$(basename $w)
                        echo $window
                        mkdir -p $window
                        cd $window
                        echo "I am here:"
                        echo $PWD
                        raxmlHPC-PTHREADS -T 6 -b 12345 -p 13579 -# 100 -m GTRCAT -s $w -n $window
                        cd ~/researchMukul/pgtr/nov8_2020/window_analysis_60perc/03_bootstrap_trees/$dir/$param/$i
                    fi
                done
                break 3
            fi
        done
    fi
done

