#/bin/bash
start=$(date +%s.%N)

python diploSHIC.py fvecSim diploid /home/khoih/fwdpy11test/data/test2.txt /home/khoih/fwdpy11test/data/test2.txt.diploid.fvec 


dur=$(echo "$(date +%s.%N) - $start" | bc)

printf "Execution time: %.6f seconds" $dur

