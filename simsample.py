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
parser.add_argument("-m", "--mutrate", help="mutrate_s", type=float)
parser.add_argument("-s", "--seed",help="seed", type=int)
parser.add_argument("-n","--nreps",help="number of replicates",type=int,default=100)
parser.add_argument("-f", "--filename", help="output_file_name")

args = parser.parse_args(sys.argv[1:])

if args.filename is None:
    raise ValueError("output file name not specified")

if args.mutrate is None:
    raise ValueError("mutation rate not specified")

if args.seed is None:
    raise ValueError("seed not specified")
	
if os.path.isfile(args.filename):
    raise Warning("output file",args.output,"exits--exiting!")
    sys.exit(0)
    
N = 1000
theta = 1000.0  # float, not int
rho = 1000.0    # float, not int


rng = fp11.GSLrng(args.seed)

for rep in range(args.nreps):
    t2f = fp11qt.GSSmo([(0,0,1),(10*N,1,1)])
    pop = fp11.SlocusPop(N)
    p = {'nregions':[fp11.Region(0,11,1.0)], #nregions must be a region
    'sregions':[fp11.GaussianS(5,6,1,0.25)],
    'recregions':[fp11.Region(0,11,1)],
    'mutrate_n':theta/float(4*N),
    'mutrate_s':args.mutrate,   # Add this as you are using mutrate_n and recrate
    'recrate':rho/float(4*N),
    # 'rates':(2.5e-4,1e-3,5e-3),  # Do not need this if using mutrate_n and recrate
    # Change to N instead of 10
    'demography':np.array([N]*(10*N+100),dtype=np.uint32),
    'gvalue':fp11tv.SlocusAdditiveTrait(2.0),
    'trait_to_fitness':t2f,
    # Do not remove selected fixations
    # from pop during simulation:
    'prune_selected': False
    }

    params = fp11.model_params.SlocusParamsQ(**p)

    fp11qt.evolve(rng,pop,params)

    # Write the pickled output:
    #The mode argument can be any of "r", "rb", "w", "wb", "x", "xb", "a" or "ab" for binary mode, or "rt", "wt", "xt", or "at" for text mode. The default is "rb".
    #"r" for reading (default), "w" for overwriting, "x" for exclusive creation, or "a" for appending

    with lzma.open(args.filename + ".lzma","ab") as f:
            pickle.dump(pop,f)
    s = pop.sample(rng=rng,nsam=2000,separate=False)
    ms = libsequence.polytable.SimData(s)
    print(str(ms), file=open(args.filename +".txt","a"))
    
