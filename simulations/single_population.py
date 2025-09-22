#!/usr/bin/env python3

#Import necessary packages
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

### Run simulation
ts = msprime.sim_ancestry(samples=1000, population_size=popsize, random_seed=seed, sequence_length=length, recombination_rate=1.8e-8, model=msprime.StandardCoalescent())

### Set filename and dump to file
simname = "/rds/general/user/jjr18/home/paper/results/single_population_ne" + str(popsize) + "_length" + str(length) + "_seed" + str(seed) + ".ts"
ts.dump(simname)
