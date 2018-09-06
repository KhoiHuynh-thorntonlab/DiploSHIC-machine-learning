#!/bin/bash

for i in $(seq 0 21);
do
	echo "0.5window$i |" >> /home/khoih/0.5total.summary.txt
	cut -f 5 /home/khoih/backup/0.5window$i.txt | sort | uniq -c >> /home/khoih/0.5window$i.summary.txt
	cat /home/khoih/0.5window$i.summary.txt >> /home/khoih/0.5total.summary.txt
done


for i in $(seq 0 99);
do
	echo "0.1window$i |" >> /home/khoih/0.1total.summary.txt
	cut -f 5 /home/khoih/backup/0.1window$i.txt | sort | uniq -c >> /home/khoih/0.1window$i.summary.txt 
	cat /home/khoih/0.1window$i.summary.txt >> /home/khoih/0.1total.summary.txt
done
