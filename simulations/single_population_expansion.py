#!/usr/bin/env python3

### Import necessary packages
import msprime
import tskit
from argparse import ArgumentParser

### Set up our command line argument parser
parser = ArgumentParser()
parser.add_argument("seed")
parser.add_argument("recent_popsize")
parser.add_argument("anc_popsize")
parser.add_argument("length")

### Set variables from the command line
args = parser.parse_args()
seed = int(args.seed)
recent_popsize = int(args.recent_popsize)
anc_popsize = int(args.anc_popsize)
length = int(args.length)


###Set up demography with expansion
# Set parameters for expansion
recent_pop = recent_popsize
anc_pop = anc_popsize
time_of_change = 10000

# Demographic model with expansion
demography = msprime.Demography()
demography.add_population(name="A", initial_size=recent_pop)

# At 10000 generations ago, population changes instantaneously
demography.add_population_parameters_change(
    time=time_of_change,
    initial_size=anc_pop,
    population="A"
)

### Sanity/debug printing of demography and simulation call
print(demography.debug())
print("ts = msprime.sim_ancestry(samples=1000, demography=demography, random_seed=seed, sequence_length=length, recombination_rate=1.8e-8, model=msprime.StandardCoalescent())")

### Run simulation
ts = msprime.sim_ancestry(samples=1000, demography=demography, random_seed=seed, sequence_length=length, recombination_rate=1.8e-8, model=msprime.StandardCoalescent())

### Set filename and dump to file
simname = "/rds/general/user/jjr18/projects/evolgenetics/live/josh/paper/results/expansion/expansion_ne" + str(recent_pop)+ "_" + str(anc_pop) + "_length" + str(length) + "_seed" + str(seed) + ".ts"
ts.dump(simname)
