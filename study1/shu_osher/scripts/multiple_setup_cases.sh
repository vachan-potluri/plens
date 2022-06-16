for N in 1 2 3 5
do
    for dof in 200 400 800
    do
        declare -i ny=$N+1
        declare -i freq=$dof/2
        dir_name=../run/dof$dof\_$ny\_$ny\_N$N\_chandrashekhar
        echo Setting up in $dir_name
        python3 setup_case.py $dir_name -f Chandrashekhar -k 4 -w $freq -r result_16Jun2022
        wait
    done
done