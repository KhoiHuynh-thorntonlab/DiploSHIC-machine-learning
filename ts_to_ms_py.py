import fwdpy11
import fwdpy11.ts
import gzip
import argparse
import sys
import numpy as np


# set number of loci and locus boundary:
NLOCI = 10
LOCUS_BOUNDARIES = [(i, i + 11) for i in range(0, NLOCI * 11, 11)]

# creat parser for command line argument for the script variable:

def make_parser():
    parser = argparse.ArgumentParser()
    # TODO: Khoi: add the help lines for these
    parser.add_argument('-filename','--filename', type=str, help = "inputfilename. Simulated population is pickled to this file")
    parser.add_argument('-theta','--theta', type=float, help="Scaled mutation rate, theta=4Nu")
    parser.add_argument('-nsam','--nsam', type=int,help="number of sampled individuals")
    parser.add_argument('-seed','--seed', type=int, help="random number seed")
    return parser

# this function to validate argument 
# input
def validate_arguments(args):
    if args.filename is None:
        raise ValueError("input file name (filename) not specified")
    if args.nsam is None:
        raise ValueError("number of sample (nsam) not specified")
    if args.theta is None:
        raise ValueError("theta value (theta) not specified")
    if args.seed is None:
        raise ValueError("seed for rng (seed) not specified")


# this is to reformat data from tree sequence sample data to a 
# MS format data 
def reformat_data(window_data, pos_window):
    
# this line is str.format() that return value in format() into 
# str {}.
# our script create a line of 1 2 3
# follows by nsetsites with position window parsing in to {}
    ms = "1 2 3\n//\nsegsites: {}\n".format(len(pos_window))

    if len(pos_window) == 0:
        return ms

    ms = ms + "positions: "
    for p in pos_window:
        ms = ms + "{} ".format(p)

    ms = ms + "\n"

    for i in range(window_data.shape[1]):
        #window_data[:,i] spliced window_data array and keep only i column
        # with seperator (delimiter) set as '' to have no space.
        # then that array is converted into string via array2string
        # which is further classified as string with str. 
        # max_line_width set numpy printing threshold for line to double wind_data shape size
        temp = str(np.array2string(window_data[:, i], separator='',max_line_width=2 * window_data.shape[0]))[1:-1]
        # NOTE: the next two lines are due to annoyances
        # with np.array2string!!!!!
        # replace any new line and space with no space via ''
        temp = temp.replace('\n', '')
        temp = temp.replace(' ', '')
        # input temp as new line after ms
        ms = ms + "{}\n".format(temp)

    return ms

# write ms format to file:
def write_ms_format(data,timepoint, pos, outfile_stub, window_size=1, step_size=0.5):
    # locus is counter and start_stop is object retrieve from locus_boundaries.
    # enumerate is to number object in locus bounderies list. Said number
    # is locus.
    for locus, start_stop in enumerate(LOCUS_BOUNDARIES):
        # np.arrange generate number from the first value from
        # start_stop[0] to last value start_stop[1] with increasing step
        # size as step_size. 
        window_starts = np.arange(start_stop[0], start_stop[1], step_size)
        for window, ws in enumerate(window_starts):
            # window_data_indexes gnerate index for all entires
            # that has position bigger than window value from wind_start
            # and smaller than window value + window_size. 
            #[0] just to get the first column of generated index.
            window_data_indexes = np.where(
                (pos >= ws) & (pos < ws + window_size))[0]
            pos_window = pos[window_data_indexes]
            # for all data with window_data_indexes, take all 
            # column
            window_data = data[window_data_indexes, :]
            # convert data into ms format using
            # reformat_data function
            ms = reformat_data(window_data, pos_window)
            print(ms)
            with open(outfile_stub+".window"+ str(window)+".timepoint" +str(timepoint)+".msdata.txt", "w") as f:
                f.write(ms)

def write_metadata(t,mean_trait,genetic_variance,mean_fitness,outfile_stub):
    with open(outfile_stub+".metadata.txt","a") as f:
        f.write(str(t)+"\t"+str(mean_trait)+"\t"+str(genetic_variance)+"\t"+str(mean_fitness)+"\n")                
                
def process_replicate(filename, repid, seed, nsam):

    with gzip.open(filename, 'rb') as f:
        pop = fwdpy11.SlocusPop.load_from_pickle_file(f)

    rng = fwdpy11.GSLrng(seed)
    np.random.seed(seed)

    # Add neutral mutations
    nm = fwdpy11.ts.infinite_sites(
        rng, pop, float(NLOCI) * args.theta / (4 * pop.N))

    # Get the node table for the pop and the
    # metadata for ancient samples
    nodes = np.array(pop.tables.nodes, copy=False)
    amd = np.array(pop.ancient_sample_metadata, copy=False)
    # These are the times for each ancient sample
    amt = nodes['time'][amd['nodes'][:, 0]]

    for t in np.unique(amt):
        # Get the indexes of the metadata corresponding
        # to this time point
        sample_indexes_at_time = np.where(amt == t)[0]

        # Calc some simple stats about the overall pop'n
        # You can figure out where to record these.
        # calculate mean genetic value for sample with
        # sample_indexes_at_time index by .mean()
        mean_trait_value = amd['g'][sample_indexes_at_time].mean()
        # calulate genetic variance for all sample by .var()
        vg = amd['g'][sample_indexes_at_time].var()
        # calculate mean fitness value w with .mean()
        wbar = amd['w'][sample_indexes_at_time].mean()
        # choose random choice:
        random_sample = np.random.choice(
            sample_indexes_at_time, nsam, replace=False)
        # flatten sample array
        samples = amd['nodes'][random_sample].flatten()
        # simplfy tree sequence array sample into tuple
        tables, idmap = fwdpy11.ts.simplify(pop, samples)
        # create data matrix from data table
        gm = fwdpy11.ts.data_matrix_from_tables(
            tables, pop.mutations, idmap[samples], True, True)
        # create array for position from neutral and selected position
        pos = np.array(gm.neutral.positions + gm.selected.positions)
        sorted_pos_indexes = np.argsort(pos)
        # input data from neutral and create array
        all_sites = np.array(gm.neutral, copy=False)
        if len(gm.selected.positions) > 0:
            # concatenate neutral and selected data
            all_sites = np.concatenate(
                (all_sites, np.array(gm.selected, copy=False)))
        # take all data according to sorted_pos_indexes
        # Re-order the matrix by the sorted positions
        all_sites = all_sites[sorted_pos_indexes, :]
        pos = pos[sorted_pos_indexes]
        write_ms_format(all_sites,t, pos, filename, 1, 0.5)
        write_metadata(t,mean_trait_value,vg,wbar,filename)

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    validate_arguments(args)
    with open(args.filename+".metadata.txt","w") as f:
        f.write("timepoint"+"\t"+"mean_trait"+"\t"+"genetic_variance"+"\t"+"mean_fitness"+"\n")
    process_replicate(args.filename, 1, args.seed, args.nsam)
