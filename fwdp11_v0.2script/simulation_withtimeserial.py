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
# They are very efficient objects in terms of speed, etc.
Key = namedtuple('Key', ['pos', 'esize', 'origin'])

DBrecord = namedtuple(
    "DBrecord", ['pos', 'esize', 'origin', 'generation', 'count'])


class Traj(object):
    """
    Track the frequency trajectory of all selected variants
    """
  
    def __init__(self, timepoints, pop, repargs):
        # We will track a dictionary of (mutation, [frequency])
        # data.  The frequency vector will be a plain list of
        # counts, which are integers.  In Python, this leads to a
        # C-like array under the hood instead of a linked list.
        # The array is MUCH better for performance.
        filename, popu, theta, rho, mutrate, seed ,repid = repargs
        self.filename = filename
        self.seed = seed
        self.popu = popu 
        self.theta = theta
        self.rho = rho
        self.mutrate = mutrate
        self.repid = repid
        self.rng=fwdpy11.GSLrng(self.seed)
        self.timepoints = timepoints
        self.data = dict()
        self.timepoints = timepoints
    def take_sample(self, pop):
            N=self.popu
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
            self.sample= samples

    def trajectory(self, pop):
            N = self.popu
            if self.timepoints < 10*pop.N and pop.generation % 100 == 0:
            # Remove any extinct variants.
            # This assumes the optimum shift happens at 10N generations.
            # We only care about variants that either "cross" the shift
            # or arise after it.  This loop removes all variants
            # arising before the shift but not crossing it.
            # We only do this every 100 generations, so that
            # we don't slow things down too much.
                nremoved = 0
                badkeys = []
                for k, v in self.data.items():
                    if k.origin + len(v) < pop.generation:
                        badkeys.append(k)
                        nremoved += 1
                for bk in badkeys:
                    self.data.pop(bk)

        # Update the mutation frequencies

        # First, get numpy arrays of the
        # population data.
            ma = np.array(pop.mutations.array())
            mc = np.array(pop.mcounts)

            for i in range(len(mc)):
            # Check that mutation is not extinct
            # and that it is NOT neutral:
                if ma['neutral'][i] == 0 and mc[i] > 0:
                # We need a unique "key" for a mutation.
                # Here, we use (position, effect size,
                # and generation when it arose).
                    k = Key(ma['pos'][i], ma['s'][i], ma['g'][i])
                    if k in self.data:
                    # if the key exists, and the mutation
                    # is not already recored as "fixed",
                    # update the frequency vector
                        if self.data[k][-1] < 2*self.popu:
                            self.data[k].append(mc[i])
                    elif mc[i] < 2*N:
                    # The mutation does not exist,
                    # so we create a frequency vector
                    # with its frequency.
                    # By definition, thus mutation
                    # must be new, meaning frequency == 1,
                    # so we assert that that is the case here:
                    # Note: the 'elif' above is sutble.
                    # We are pruning fixations that do not
                    # cross the optimum shift, but they are
                    # still in the population, so we need
                    # to avoid re-adding them over and over.
                        assert mc[i] == 1, "mutation frequency error {} {} {}".format(
                            mc[i], ma['neutral'][i], ma['g'][i])
                    # Insert the new record:
                        self.data[k] = [mc[i]]

    def __call__(self, pop):

        N = self.popu
        self.timepoints = pop.generation
        self.trajectory(pop)
        if self.timepoints >10*pop.N:
                    with lzma.open(self.filename +'.'+str(self.repid)+'.'+ str(self.timepoints) + '.lzma', "ab") as f:
                        print("pickling")
                        pickle.dump((pop), f, -1)
                        print("done pickling")
                    print("writing samples")

#demograph generation total is 10*popu+100
def runsim(repargs):
    filename, popu, theta, rho, mutrate, seed ,repid = repargs
    pop = fp11.SlocusPop(popu)
    rng = fp11.GSLrng(seed)
    optima = fwdpy11.genetic_values.GSSmo([(0, 0, 1), (10*popu, 1, 1)])
    gv = fwdpy11.genetic_values.SlocusAdditive(2.0, optima)
    p = {'nregions': [fp11.Region(0,11,1.0)],
         'sregions': [fp11.GaussianS(5, 6, 1, 0.25)],
         'recregions': [fp11.Region(0, 11, 1)],
         # No neutral mutations in this example
         'rates':(theta/float(4*popu), mutrate, rho/float(4*popu)),
         'gvalue': gv,
         'prune_selected': False,
         'demography': np.array([popu]*(10*popu+100), dtype=np.uint32)
         }
    timepoints = pop.generation
    sampler = Traj(timepoints, pop, repargs)
    params = fp11.model_params.ModelParams(**p)
    fwdpy11.wright_fisher.evolve(rng, pop, params, sampler)
    return repid, pop, sampler.data 



def make_trajectory_DataFrame(data, repid, popu):
    """ 
    Take the frequency data from our sampler
    and coerce it into a DataFrame object
    that we can then write out to an sqlite3
    database.
    We want to end up with a data frame with 
    the following columns:
    repid generation pos esize origin count
    repid = simulation replicate
    generation = generation in the simulation
    pos = mutation position
    esize = effect size
    origin = generation when mutation first arose
    count = number of time the mutation is present in this generation.
    Thus, the DataFrame will completely specify the entire mutation
    frequency history for all mutations in the replicate.
    """
    dfdata = []
    for k, v in data.items():
        if k.origin + len(v) - 1 >= 10*popu:
            # Only keep mutations
            # whose trajectories end after the optimum shift.
            t = [DBrecord(k.pos, k.esize, k.origin, k.origin+j, v[j])
                 for j in range(len(v))]
            dfdata.extend(t)
    # The nice thing about namedtuple object is
    # that a list of them can be the data
    # for a DataFrame.  The "trick" is to
    # use the _fields variable as the column names!
    # These field names are what is defined above on line 18.
    df = pd.DataFrame(dfdata, columns=DBrecord._fields)
    # Don't forget to add the replicate id to
    # the dataframe!!!!
    df['repid'] = [repid]*len(df.index)

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mutrate", help="mutrate_s", type=float)
    parser.add_argument("-s", "--seed", help="seed", type=int)
    parser.add_argument("-n", "--nreps", help="number of replicates", type=int, default=100)
    parser.add_argument("-f", "--filename", help="output_file_name")
    parser.add_argument("-t", "--theta",help="theta value", type=int)
    parser.add_argument("-r", "--rho",help="rho value", type=int)
    parser.add_argument("-p", "--popu",help="population value", type=int)
    args = parser.parse_args(sys.argv[1:])
    

    if args.filename is None:
        raise ValueError("output file name not specified")

    if args.mutrate is None:
        raise ValueError("mutation rate not specified")

    if args.seed is None:
        raise ValueError("seed not specified")
    if args.popu is None:
        raise ValueError("population value not specified")
    if args.theta is None:
        raise ValueError("theta value not specified")
    if args.rho is None:
        raise ValueError("rho value not specified")    

    if args.nreps is None:
        raise ValueError("number of rep not specified")    

    if os.path.isfile(args.filename):
        raise Warning("output file", args.output, "exits--exiting!")
        sys.exit(0)

    np.random.seed(args.seed)

    # generate seeds for each replicate
    repseeds = np.random.choice(range(1000000), args.nreps, replace=False)

    # create a list of argument tuples for each replicate
    repargs = [((args.filename), (args.popu), (args.theta), (args.rho), (args.mutrate), (i), (j))
               for i, j in zip(repseeds, range(args.nreps))]


    # delete the pickled pop file if it exists,
    # which is important b/c we are going to
    # use this file in append mode.
    if os.path.exists(args.filename + '.lzma'):
        os.remove(args.filename + '.lzma')
    if os.path.exists(args.filename + '.data.txt'):
        os.remove(args.filename + '.data.txt')
    with concurrent.futures.ProcessPoolExecutor() as pool:
        futures = {pool.submit(runsim, i) for i in repargs}
        for fut in concurrent.futures.as_completed(futures):
            result = fut.result()
            repid, pop, Traj.trajectory.data = result
            # The output file name here is hard-coded.
            # This is bad, and should be an option.
            with sqlite3.connect(args.filename + ".sqlite3") as conn:
                df = make_trajectory_DataFrame(Traj.trajectory.data, repid, args.popu)
                df.to_sql('data', conn, if_exists='append',
                          index=False, chunksize=5000)

