#!/usr/bin/env python3

### Import necessary packages
import msprime
import tskit
from argparse import ArgumentParser

### Set up our command line argument parser
parser = ArgumentParser()
parser.add_argument("seed")
parser.add_argument("popsize")
parser.add_argument("length")

### Set variables from the command line
args = parser.parse_args()
seed = int(args.seed)
popsize = int(args.popsize)
length = int(args.length)


###Set up demography with island model
demography = msprime.Demography.island_model(initial_size=[popsize, popsize, popsize, popsize, popsize], migration_rate=0.01)


### Sanity/debug printing of demography
print(demography.debug())

### Set up samples to take from each deme
samples = {"pop_0": 200, "pop_1": 200, "pop_2": 200, "pop_3": 200, "pop_4": 200,}

### Run simulation
ts = msprime.sim_ancestry(samples=samples, demography=demography, random_seed=seed, sequence_length=length, recombination_rate=1.8e-8, model=msprime.StandardCoalescent())

### Set filename and dump to file
simname = "/rds/general/user/jjr18/projects/evolgenetics/live/josh/paper/results/island/island_ne" + str(anc_pop) + "_length" + str(length) + "_seed" + str(seed) + ".ts"
ts.dump(simname)
