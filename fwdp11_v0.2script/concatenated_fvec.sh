#!/bin/bash

list=( 10 )

for i in "${list[@]}";
do
       for j in $(seq 10001 10100);
        do
                for k in $(seq 0 109);
                do
			paste /home/khoih/stat.txt highmutrate/high.${j}.rep${i}.window${k}.txt.haploid.fvec > highmutrate/withcol/high.${j}.rep${i}.window${k}.haploid.fvec
			tail -n +2 highmutrate/withcol/high.${j}.rep${i}.window${k}.haploid.fvec >> temporary/highmutratestat.$i.$j.txt
                done
        done
done

