for w in $1/window*; do
    ls -lh $w
    if [[ $w == *.reduced ]]
    then
	echo 'pass'
    else
	window=$(basename $w)
	echo $window
	raxmlHPC-PTHREADS -T 6 -b 12345 -p 13579 -# 100 -m $2 -s $w -n $window
	split -l 1 --additional-suffix=."$window"_bs RAxML_bootstrap.$window
    fi
done

