for w in $1/window*; do
    ls -lh $w
    window=$(basename $w)
    echo $window
    raxmlHPC-PTHREADS -T 6 -f a -x 12345 -p 12345 -s $w -n $window -m $2 -#100 
done

