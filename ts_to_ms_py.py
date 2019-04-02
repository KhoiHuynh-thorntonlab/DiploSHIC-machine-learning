import fwdpy11
import fwdpy11.ts
import gzip
import argparse
import sys
import numpy as np

NLOCI = 10
LOCUS_BOUNDARIES = [(i, i + 11) for i in range(0, NLOCI * 11, 11)]


def make_parser():
    parser = argparse.ArgumentParser()

    # TODO: Khoi: add the help lines for these
    parser.add_argument('-filename','--filename', type=str, help = "inputfilename. Simulated population is pickled to this file")
    # I keep theta default at 1000 as it is the default we normally used.
    parser.add_argument('-theta','--theta', type=float, help="Scaled mutation rate, theta=4Nu")
    parser.add_argument('-nsam','--nsam', type=int,help="number of sampled individuals")
    parser.add_argument('-seed','--seed', type=int, help="random number seed")

    return parser


def validate_arguments(args):
    if args.filename is None:
        raise ValueError("input file name (filename) not specified")
    if args.nsam is None:
        raise ValueError("number of sample (nsam) not specified")
    if args.theta is None:
        raise ValueError("theta value (theta) not specified")
    if args.seed is None:
        raise ValueError("seed for rng (seed) not specified")


def reformat_data(window_data, pos_window):
    ms = "1 2 3\n//\nsegsites: {}\n".format(len(pos_window))

    if len(pos_window) == 0:
        return ms

    ms = ms + "positions: "
    for p in pos_window:
        ms = ms + "{} ".format(p)

    ms = ms + "\n"

    for i in range(window_data.shape[1]):
        temp = str(np.array2string(window_data[:, i], separator=''))[1:-1]
        # NOTE: the next two lines are due to annoyances
        # with np.array2string!!!!!
        temp = temp.replace('\n', '')
        temp = temp.replace(' ', '')
        ms = ms + "{}\n".format(temp)

    return ms


def write_ms_format(data, timepoint, pos, outfile_stub, window_size=1, step_size=0.5):
    for locus, start_stop in enumerate(LOCUS_BOUNDARIES):
        window_starts = np.arange(start_stop[0], start_stop[1], step_size)
        for window, ws in enumerate(window_starts):
            window_data_indexes = np.where(
                (pos >= ws) & (pos < ws + window_size))[0]
            pos_window = pos[window_data_indexes]
            window_data = data[window_data_indexes, :]
            ms = reformat_data(window_data, pos_window)
            np.set_printoptions(threshold=sys.maxsize)
            print(ms)
            with open(outfile_stub+".window"+ str(window)+".timepoint" +str(timepoint)+".msdata.txt", "w") as f:
                f.write(ms)
def write_metadata(t,mean_trait,genetic_value,mean_fitness,outfile_stub):
    with open(outfile_stub+".metadata.txt","a") as f:
        f.write(str(t)+"\t"+str(mean_trait)+"\t"+str(genetic_value)+"\t"+str(mean_fitness)+"\n")

def process_replicate(filename, repid, seed, nsam):

### Open simulated population :
    with gzip.open(filename +".gz", 'rb') as f:
        pop = fwdpy11.SlocusPop.load_from_pickle_file(f)

### RNG seed generation

    rng = fwdpy11.GSLrng(seed)
    np.random.seed(seed)

# Add neutral mutations
# since tree sequence simulation do not
# involve neutral mutation. it only
# simulate selected mutation.

    nm = fwdpy11.ts.infinite_sites(
        rng, pop, float(NLOCI) * args.theta / (4 * pop.N))

# Get the node table for the pop and the
# metadata for ancient samples
    nodes = np.array(pop.tables.nodes, copy=False)
    amd = np.array(pop.ancient_sample_metadata, copy=False)
# These are the times for each ancient sample
# amt is time element from ancient sample metadata in node table:
    amt = nodes['time'][amd['nodes'][:, 0]]
# t is the timepoint of the sample
    for t in np.unique(amt):
        # Get the indexes of the metadata corresponding
        # to this time point
        # np.where rereturn indexes of entries that satisfy amt == t
        sample_indexes_at_time = np.where(amt == t)[0]

        # Calc some simple stats about the overall pop'n
        # You can figure out where to record these.
        # retrieve mean trait value g from ancient metadata from node table with index
        # specified by amt (timepoint)
        mean_trait_value = amd['g'][sample_indexes_at_time].mean()
        vg = amd['g'][sample_indexes_at_time].var()
        wbar = amd['w'][sample_indexes_at_time].mean()

        # randomly choose entries that has index classfied by sample_indexes_at_time and condition amt==t
        random_sample = np.random.choice(
            sample_indexes_at_time, nsam, replace=False)
        samples = amd['nodes'][random_sample].flatten()
        tables, idmap = fwdpy11.ts.simplify(pop, samples)

        gm = fwdpy11.ts.data_matrix_from_tables(
            tables, pop.mutations, idmap[samples], True, True)
        pos = np.array(gm.neutral.positions + gm.selected.positions)
        sorted_pos_indexes = np.argsort(pos)
        all_sites = np.array(gm.neutral, copy=False)
        if len(gm.selected.positions) > 0:
            all_sites = np.concatenate(
                (all_sites, np.array(gm.selected, copy=False)))
        # Re-order the matrix by the sorted positions
        all_sites = all_sites[sorted_pos_indexes, :]
        pos = pos[sorted_pos_indexes]
        # convert sample to ms format with write_ms_format function
        write_ms_format(all_sites,t, pos, filename, 1, 0.5)
        write_metadata(t,mean_trait_value,vg,wbar,filename)

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    validate_arguments(args)
    with open(args.filename+".metadata.txt","w") as f:
        f.write("timepoint"+"\t"+"mean_trait"+"\t"+"genetic_value"+"\t"+"mean_fitness"+"\n")
    process_replicate(args.filename, 1, args.seed, args.nsam)

