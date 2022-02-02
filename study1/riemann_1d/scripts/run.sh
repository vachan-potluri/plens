# A bash script to run simulations for all values of N for a specific case

for N in 1 2 3 5
do
    cd dof400_12_12_N$N\_hllc
    # echo $N
    nohup mpirun -np 4 ../../../../../build/plens.out 1 $N > result/plens.log 2>&1 &
    wait
    cd ..
done
