#!/usr/bin/env python3

### Import packages
import msprime
import numpy as np
from IPython.display import SVG
import tskit
import pandas as pd
import re
from argparse import ArgumentParser

### Set up our command line argument parser
parser = ArgumentParser()
parser.add_argument("sim")
parser.add_argument("type")

### Set variables from the command line
args = parser.parse_args()
sim = args.sim
type = args.type

### Define path to tree sequence - change as needed -  and then load the sequence
ts_path = "/rds/general/user/jjr18/results/" + str(type) + "/" + str(sim) + ".ts_mutated.ts"
ts_mutated = tskit.load(ts_path)

### Identify recurrent mutations
# Set up lists to hold positions of recurrent mutations
same_base_sites = []
different_base_sites = []

# Iterate through tree sequence
tree_iter = ts_mutated.trees()
tree = next(tree_iter)
tree.seek(0)

for site in ts_mutated.sites():
    # Advance tree if needed
    while site.position >= tree.interval.right:
        tree = next(tree_iter)

    if len(site.mutations) == 2:
        derived_states = []

        for mut in site.mutations:
            count = sum(1 for _ in tree.samples(mut.node))
            if count == 1:
                derived_states.append(mut.derived_state)
        # Determine if both mutations are to the same nucleotide (i.e. a doubleton) or not
        if len(derived_states) == 2:
            pos = int(site.position)
            if derived_states[0] == derived_states[1]:
                same_base_sites.append(pos)
            else:
                different_base_sites.append(pos)


### Write same base sites (a.k.a doubletons) to file
filename = "/rds/general/user/jjr18/results/" + str(type) + "/recurrent_sites_"+ str(sim) +".txt"
with open(filename, "w") as f:
    f.write(",".join(map(str, same_base_sites)))
