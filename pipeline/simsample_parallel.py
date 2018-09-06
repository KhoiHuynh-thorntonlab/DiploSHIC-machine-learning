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
import libsequence.polytable
# KRT: need the windows module to
# split up data
import libsequence.windows
import concurrent.futures


# Python style guidelines
# recommend ALLCAPS for global
# variables
N = 1000
THETA = 1100.0
RHO = 1100.0


def take_sample(rng, pop):
    s = pop.sample(rng=rng, nsam=25, separate=False)

    ms = libsequence.polytable.SimData(s)
    w = libsequence.windows.Windows(ms, window_size=1.0, step_len=0.1,
                                    starting_pos=0.0, ending_pos=11.0)

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
    return samples


# Functions run via concurrent futures are simplest to
# implement if they take a tuple of parameters
def runsim(args):
    mutrate, seed, repid = args
    rng = fp11.GSLrng(seed)
    t2f = fp11qt.GSSmo([(0, 0, 1), (10*N, 1, 1)])
    pop = fp11.SlocusPop(N)
    p = {'nregions': [fp11.Region(0, 11, 1.0)],  # nregions must be a region
         'sregions': [fp11.GaussianS(5, 6, 1, 0.25)],
         'recregions': [fp11.Region(0, 11, 1)],
         'mutrate_n': THETA/float(4*N),
         'mutrate_s': mutrate,
         'recrate': RHO/float(4*N),
         'demography': np.array([N]*(10*N+100), dtype=np.uint32),
         'gvalue': fp11tv.SlocusAdditiveTrait(2.0),
         'trait_to_fitness': t2f,
         'prune_selected': False
         }

    params = fp11.model_params.SlocusParamsQ(**p)

    fp11qt.evolve(rng, pop, params)
    samples = take_sample(rng, pop)
    return repid, pop, samples


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mutrate", help="mutrate_s", type=float)
    parser.add_argument("-s", "--seed", help="seed", type=int)
    parser.add_argument(
        "-n", "--nreps", help="number of replicates", type=int, default=100)
    parser.add_argument("-f", "--filename", help="output_file_name")

    args = parser.parse_args(sys.argv[1:])

    if args.filename is None:
        raise ValueError("output file name not specified")

    if args.mutrate is None:
        raise ValueError("mutation rate not specified")

    if args.seed is None:
        raise ValueError("seed not specified")

    if os.path.isfile(args.filename):
        raise Warning("output file", args.output, "exits--exiting!")
        sys.exit(0)

    np.random.seed(args.seed)

    # generate seeds for each replicate
    repseeds = np.random.choice(range(1000000), args.nreps, replace=False)

    # create a list of argument tuples for each replicate
    repargs = [(args.mutrate, i, j)
               for i, j in zip(repseeds, range(args.nreps))]

    # delete the pickled pop file if it exists,
    # which is important b/c we are going to 
    # use this file in append mode.
    if os.path.exists(args.filename + '.lzma'):
        os.remove(args.filename + '.lzma')

    with concurrent.futures.ProcessPoolExecutor() as pool:
        futures = {pool.submit(runsim, i) for i in repargs}
        for fut in concurrent.futures.as_completed(futures):
            result = fut.result()
            repid, pop, samples = result
            print("{} done, {} {}".format(repid, pop.generation, len(samples)))
            with lzma.open(args.filename + '.lzma', "ab") as f:
                print("pickling")
                pickle.dump((repid, pop), f, -1)
                print("done pickling")
            print("writing samples")
            for i, j in zip(samples, range(len(samples))):
                with open(args.filename + ".rep{}.window{}.txt".format(repid, j), "w") as f:
                    # Write a fake discoal header line with fake random number seeds
                    f.write("discoal 50 1 1000 -t 1000 -r 1000\n1 2 3\n")
                    f.write(str(i[1]))
            print("done writing samples")
