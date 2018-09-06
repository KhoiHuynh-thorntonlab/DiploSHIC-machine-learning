#/bin/bash

start=$(date +%s.%N)


for i in $(seq 0 99)

do
	for j in $(seq 0 99)
	do

		python diploSHIC.py fvecSim haploid /home/khoih/test/trial.rep$i.window$j.txt /home/khoih/test/stat/trial.rep$i.window$j.haploid.fvec --numSubWins 11 
	done
done


for i in $(seq 0 99);
do
	for j in $(seq 0 99);
	do
		paste /home/khoih/stat.txt /home/khoih/test/stat/trial.rep$i.window$j.haploid.fvec > /home/khoih/test/withcol/trial.rep$i.window$j.haploid.fvec
	done
done


for i in $(seq 0 99);
do
	for j in $(seq 0 99);
	do
		python /home/khoih/diploSHIC/diploSHIC.py predict /home/khoih/dis/modelfile.txt.json /home/khoih/dis/modelfile.txt.weights.hdf5 /home/khoih/test/withcol/trial.rep$i.window$j.haploid.fvec /home/khoih/test/result/trial.rep$i.window$j.result.txt
	done
done


dur=$(echo "$(date +%s.%N) - $start" | bc)

printf "Execution time: %.6f seconds" $dur

