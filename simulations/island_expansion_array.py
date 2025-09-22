#!/usr/bin/env python3

### Import necessary packages
import msprime
import numpy as np
import tskit
import pandas as pd
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
popsize = int(args.popsize)
length = int(args.length)

###Set up demography with island model and expansion
# Set parameters for expansion
recent_pop = recent_popsize
anc_pop = anc_popsize
time_of_change = 10000

#Island model with 5 demes with expansion
demography = msprime.Demography.island_model(initial_size=[recent_pop, recent_pop, recent_pop, recent_pop, recent_pop], migration_rate=0.01)

# At 10000 generations ago, population stops changing
demography.add_population_parameters_change(
    time=time_of_change,
    initial_size=anc_pop
)

### Sanity/debug printing of demography and simulation call
print(demography.debug())
print("ts = msprime.sim_ancestry(samples=1000, demography=demography, random_seed=seed, sequence_length=length, recombination_rate=1.8e-8, model=msprime.StandardCoales")

### Set up samples to take from each deme
samples = {"pop_0": 200, "pop_1": 200, "pop_2": 200, "pop_3": 200, "pop_4": 200,}

### Run simulation
ts = msprime.sim_ancestry(samples=samples, demography=demography, random_seed=seed, sequence_length=length, recombination_rate=1.8e-8, model=msprime.StandardCoalescent())

### Set filename and dump to file
simname = "/rds/general/user/jjr18/projects/evolgenetics/live/josh/paper/results/island_expansion/island_expansion_ne" + str(recent_pop)+ "_" + str(anc_pop) + "_length" + str(length) + "_seed" + str(seed) + ".ts"
ts.dump(simname)
