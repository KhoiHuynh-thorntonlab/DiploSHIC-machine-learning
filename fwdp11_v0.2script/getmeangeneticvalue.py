import fwdpy11 as fp11
import fwdpy11.model_params
import fwdpy11.genetic_values
import fwdpy11.genetic_value_noise
import fwdpy11.wright_fisher
import numpy as np
from collections import namedtuple
import sqlite3
import pandas as pd
import concurrent.futures
import sys
import numpy as np
import pickle  # to pickle the output
import lzma    # to compress the pickled file (see fwdpy11 manual0)
import argparse
import os
import libsequence
import libsequence.polytable
# KRT: need the windows module to
# split up data
import libsequence.windows
import subprocess
# namedtuples are very convenient ways to record data.
import subprocess
import tensorflow as tf
from keras import backend as K
os.environ["CUDA_VISIBLE_DEVICES"]="0"




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename", help="output_file_name")
    parser.add_argument("-bt", "--begintime",help="first timepoints", type=int)
    parser.add_argument("-et", "--endtime",help="last timepoints", type=int)
    parser.add_argument("-s", "--seed",help="seed", type=int)
    parser.add_argument("-br", "--beginrep",help="first timepoints", type=int)
    parser.add_argument("-er", "--endrep",help="last timepoints", type=int)
    parser.add_argument("-p", "--popu",help="population value", type=int)
    args = parser.parse_args(sys.argv[1:])
    

    if args.filename is None:
        raise ValueError("output file name not specified")
    if args.beginrep is None:
        raise ValueError("begin rep value not specified")
    if args.endrep is None:
        raise ValueError("end rep value not specified")

    if args.seed is None:
        raise ValueError("seed not specified")
    if args.begintime is None:
        raise ValueError("begin timepoints value not specified")
    if args.endtime is None:
        raise ValueError("end timepoints value not specified")

    if os.path.isfile(args.filename):
        raise Warning("output file", args.output, "exits--exiting!")
        sys.exit(0)

with tf.device('/device:GPU:0'):
    for k in range(args.beginrep,args.endrep+1):
        for l in range(args.begintime,args.endtime+1):
            with lzma.open(args.filename +"."+str(k)+"."+str(l)+".lzma", "rb") as f:
                pop = pickle.load(f)
            dip_metadata = np.array(pop.diploid_metadata)
            gen_val = dip_metadata['g'].mean()
            with open(args.filename +"."+str(k)+"."+str(l)+".genval.txt", "w") as f:
                    # Write a fake discoal header line with fake random number seeds
                    f.write(str(k)+' ')
                    f.write(str(l)+' ')
                    f.write(str(gen_val))
