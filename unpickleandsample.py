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

for x in range (10):

	seed = random.randint(1,200000)
	rng = fp11.GSLrng(seed)
	with lzma.open(args.input +".lzma", "rb") as f:
		pop = pickle.load(f)
	s = pop.sample(rng=rng,nsam=100,separate=False)
	ms = libsequence.polytable.SimData(s)
	# (have to add this line so that diploshic can generate vector from the samples) print("fwdpy11 100 10", file=open(args.filename +".txt","a"))
	print(str(ms), file=open(args.filename +".txt","a")) 
      
