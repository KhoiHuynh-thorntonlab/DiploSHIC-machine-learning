import fwdpy11 as fp11
import fwdpy11.trait_values as fp11tv
import fwdpy11.wright_fisher_qtrait as fp11qt
import fwdpy11.model_params
import numpy as np
import pickle  # to pickle the output
import lzma    # to compress the pickled file (see fwdpy11 manual0)
import sys
import argparse
import os
import random
import libsequence.polytable

parser = argparse.ArgumentParser()
#parser.add_argument("-s", "--seed", help="seed value")
parser.add_argument("-i", "--input", help="file inputname")
parser.add_argument("-o", "--filename", help="file outputname")
args = parser.parse_args()


with lzma.open(args.input +".lzma", "rb") as f:
        pop = pickle.load(f)
        # (have to add this line so that diploshic can generate vector from the samples)
        #print("fwdpy11 100 10", file=open(args.filename +".txt","a"))
        #print(seed + seed, file=open(args.filename+".txt","a"))
        #print("", file=open(args.filename +".txt","a"))
print(pop.generation)
for i,j in zip(pop.fixations, pop.fixation_times):
     if i.neutral is False: # Only process selected mutatios
        if i.g < 3665:
           # This is from standing variation.
          print("S",i.pos,i.s,i.g,j)
        else:
           # A new mutation
          print ("N",i.pos,i.s,i.g,j)
