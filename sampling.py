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

# namedtuples are very convenient ways to record data.


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

    for k in range(args.beginrep,args.endrep+1):
        for l in range(args.begintime,args.endtime+1):
            with lzma.open(args.filename +"."+str(k)+"."+str(l)+".lzma", "rb") as f:
                pop = pickle.load(f)
            N=args.popu
            full_list=np.arange(0,pop.N, dtype=np.uint32)
            # the second argument in np.random.choice is the number of individual sampled:
            individuals=np.random.choice(full_list, 100)
            s = pop.sample(individuals)
    # Convert samples from version 0.2 (neutral, selected) format to (position, genotype)
    # Then sort the combined samples by position:
            neutral, selected = fwdpy11.sampling.matrix_to_sample(s)
            combined = sorted(neutral + selected, key=lambda x: x[0])
            ms = libsequence.polytable.SimData(combined)
            w = libsequence.windows.Windows(ms, window_size=1.0, step_len=0.1,starting_pos=0.0, ending_pos=11.0)

            left_end = 0.0
            samples = []
            for i in range(len(w)):
                wi = w[i]
                newpos = [wi.position(i)-left_end for i in range(wi.numsites())]
                new_window = libsequence.polytable.SimData(newpos,
                                                       [wi[i] for i in range(len(wi))])
                samples.append((left_end, new_window))
        # need to update left_end
        # to get positions correct
                left_end += 0.1
            samples = samples
            for i, j in zip(samples, range(len(samples))):
                        with open(args.filename + "." + str(l) + ".rep{}.window{}.txt".format(k, j), "w") as f:
                    # Write a fake discoal header line with fake random number seeds
                            f.write("discoal 50 1 1000 -t 1000 -r 1000\n1 2 3\n")
                            f.write(str(i[1]))
    print("done writing samples")
                 
