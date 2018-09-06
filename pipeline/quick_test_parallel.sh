#/bin/bash


# Neutral training data
mspms 100 2000 -t 100 -r 100 10000 > /home/khoih/dis/neutral_data_msprime.txt


# You have to make and maintain all the directories yourself:
if [ ! -d /home/khoih/dis/hardTraining ]
then 
    mkdir /home/khoih/dis/hardTraining
else
    rm -f /home/khoih/dis/hardTraining/*
fi

if [ ! -d /home/khoih/dis/softTraining ]
then 
    mkdir /home/khoih/dis/softTraining
else
    rm -f /home/khoih/dis/softTraining/*
fi


WIN=0
if [ -e /home/khoih/dis/simcommands.txt ]
then
    rm -f /home/khoih/dis/simcommands.txt
fi
if [ -e /home/khoih/dis/fveccommands.txt ]
then
    rm -f /home/khoih/dis/fveccommands.txt
fi
#window 0 to 4 for sweep position

for x in $(seq 0.05 0.1 0.45)
do
    # Make selected data
    echo "/home/khoih/discoal/discoal 100 2000 1000 -ws 0 -Pa 250 2500 -x $x -t 100 -r 100 > /home/khoih/dis/hard_$WIN.discoal" >> /home/khoih/dis/simcommands.txt
    echo "/home/khoih/discoal/discoal 100 2000 1000 -ws 0 -Pa 250 2500 -x $x -Pf 0.1 0.5 -t 100 -r 100 > /home/khoih/dis/soft_$WIN.discoal" >> /home/khoih/dis/simcommands.txt
    #The file names cannot have multiple _ in them
    echo "python3 /home/khoih/diploSHIC/diploSHIC.py fvecSim haploid /home/khoih/dis/hard_$WIN.discoal /home/khoih/dis/hardTraining/hard_$WIN.txt" >> /home/khoih/dis/fveccommands.txt
    echo "python3 /home/khoih/diploSHIC/diploSHIC.py fvecSim haploid /home/khoih/dis/soft_$WIN.discoal /home/khoih/dis/softTraining/soft_$WIN.txt" >> /home/khoih/dis/fveccommands.txt
    WIN=$(($WIN+1))
done

#generate 0.5 hard/soft sweep position
echo "/home/khoih/discoal/discoal 100 2000 1000 -ws 0 -Pa 250 2500 -x 0.5 -t 100 -r 100 > /home/khoih/dis/hard_5.discoal" >> /home/khoih/dis/simcommands.txt
echo "/home/khoih/discoal/discoal 100 2000 1000 -ws 0 -Pa 250 2500 -x 0.5 -Pf 0.1 0.5 -t 100 -r 100 > /home/khoih/dis/soft_5.discoal" >> /home/khoih/dis/simcommands.txt
#The file names cannot have multiple _ in them
echo "python3 /home/khoih/diploSHIC/diploSHIC.py fvecSim haploid /home/khoih/dis/hard_5.discoal /home/khoih/dis/hardTraining/hard_5.txt" >> /home/khoih/dis/fveccommands.txt
echo "python3 /home/khoih/diploSHIC/diploSHIC.py fvecSim haploid /home/khoih/dis/soft_5.discoal /home/khoih/dis/softTraining/soft_5.txt" >> /home/khoih/dis/fveccommands.txt

# window 6 to 10 for sweep position:

WIN2=6
for x in $(seq 0.55 0.1 0.95)
do
    # Make selected data
    echo "/home/khoih/discoal/discoal 100 2000 1000 -ws 0 -Pa 250 2500 -x $x -t 100 -r 100 > /home/khoih/dis/hard_$WIN2.discoal" >> /home/khoih/dis/simcommands.txt
    echo "/home/khoih/discoal/discoal 100 2000 1000 -ws 0 -Pa 250 2500 -x $x -Pf 0.1 0.5 -t 100 -r 100 > /home/khoih/dis/soft_$WIN2.discoal" >> /home/khoih/dis/simcommands.txt
    #The file names cannot have multiple _ in them
    echo "python3 /home/khoih/diploSHIC/diploSHIC.py fvecSim haploid /home/khoih/dis/hard_$WIN2.discoal /home/khoih/dis/hardTraining/hard_$WIN2.txt" >> /home/khoih/dis/fveccommands.txt
    echo "python3 /home/khoih/diploSHIC/diploSHIC.py fvecSim haploid /home/khoih/dis/soft_$WIN2.discoal /home/khoih/dis/softTraining/soft_$WIN2.txt" >> /home/khoih/dis/fveccommands.txt
    WIN2=$(($WIN2+1))
done



# Make the features from neutral data
echo "python3 /home/khoih/diploSHIC/diploSHIC.py fvecSim haploid /home/khoih/dis/neutral_data_msprime.txt /home/khoih/dis/neutral_data_msprime_features.txt" >> /home/khoih/dis/fveccommands.txt

parallel --jobs 40 < /home/khoih/dis/simcommands.txt
parallel --jobs 40  < /home/khoih/dis/fveccommands.txt


if [ ! -d /home/khoih/dis/training_out ]
then
    mkdir /home/khoih/dis/training_out
else 
    rm -f /home/khoih/dis/training_out/*
fi

python3 /home/khoih/diploSHIC/diploSHIC.py makeTrainingSets /home/khoih/dis/neutral_data_msprime_features.txt /home/khoih/dis/softTraining/soft /home/khoih/dis/hardTraining/hard 5 0,1,2,3,4,6,7,8,9,10 /home/khoih/dis/training_out
    
python3 /home/khoih/diploSHIC/diploSHIC.py train /home/khoih/dis/training_out/ /home/khoih/dis/training_out/ /home/khoih/dis/modelfile.txt

# As of Aug 22, things work for KRT up until now.

# Attempting to get a real data set analyzed is failing for me.
# My attempt is:
