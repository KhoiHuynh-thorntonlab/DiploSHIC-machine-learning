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
    parser.add_argument('filename', type=str)
    parser.add_argument('theta', type=float)
    parser.add_argument('nsam', type=int)
    parser.add_argument('seed', type=int)
    return parser

# this is to reformat data from tree sequence sample data to a 
# MS format data 
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

# write ms format to file:
def write_ms_format(data, pos, outfile_stub, window_size=1, step_size=0.5):
    for locus, start_stop in enumerate(LOCUS_BOUNDARIES):
        window_starts = np.arange(start_stop[0], start_stop[1], step_size)
        for window, ws in enumerate(window_starts):
            window_data_indexes = np.where(
                (pos >= ws) & (pos < ws + window_size))[0]
            pos_window = pos[window_data_indexes]
            window_data = data[window_data_indexes, :]
            ms = reformat_data(window_data, pos_window)
            print(ms)
            with open("foo", "w") as f:
                f.write(ms)


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
        mean_trait_value = amd['g'][sample_indexes_at_time].mean()
        vg = amd['g'][sample_indexes_at_time].var()
        wbar = amd['w'][sample_indexes_at_time].mean()

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
        write_ms_format(all_sites, pos, "foo", 1, 0.5)


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])

    process_replicate(args.filename, 1, args.seed, args.nsam)
