'''Main script to run to execute REMC simulation
'''

# Import of necessary modules

import numpy as np
import math as m
import os
import sys
import argparse

import copy
import random
import time
import multiprocessing as mp
from tqdm import tqdm
import matplotlib.pyplot as plt

from components import *

# Main program functions

def remc_full_sim(sequence,
                  protein_name, 
                  optimal_E, 
                  min_t = 160, 
                  max_t = 220, 
                  mc_steps = 500, 
                  replica_number = 25,
                  cutoff_runtime = 120,
                  nb_of_runs = 5):
    '''
    Run a full REMC simulation on a HP sequence.

    Args
    ---
    sequence 
    protein_name
    optimal_E 
    min_t 
    max_t 
    mc_steps  
    replica_number
    cutoff_runtime
    nb_of_runs

    Returns
    ---
    calculated_conformation
    '''
    # Start messages

    print(f"Running REMC simulation with the following parameters:")
    print(f"Sequence: {sequence}")
    print(f"Protein Name: {protein_name}")
    print(f"Optimal Energy: {optimal_E}")
    print(f"Min Temperature: {min_t}")
    print(f"Max Temperature: {max_t}")
    print(f"Monte Carlo Steps: {mc_steps}")
    print(f"Number of Replicas: {replica_number}")
    print(f"Cutoff Runtime: {cutoff_runtime} seconds")
    print(f"Number of Runs: {nb_of_runs}")

    # Number elements of the input sequence
    n_sequence = [f'{aa.upper()}{i + 1}' for i, aa in enumerate(sequence)]

    # Create lists of temperatures and conformations
    conf, temp = create_protein_conformations(protein_name, n_sequence, min_t, max_t, replica_number)

    # Init multiprocessing.Pool()
    print("Number of processors: ", mp.cpu_count())
    pool = mp.Pool(mp.cpu_count())

    # Create as many model instances as runs
    models = []
    for run in range(nb_of_runs):
        model = REMC(conf, temp, mc_steps)
        models.append(model)
    
    # Run the the remc sims in parallel and save final conformations for each run
    final = [pool.apply(model.run_remc_sim, args=(optimal_E, cutoff_runtime)) for model in models]

    # Close pooling
    pool.close()    
    
    # Plot and save every final plot
    for i, c in enumerate(final) :
        plot_final_conformation(c, '../results', f'Run{i}')


def main():
    '''Argument parser for the main program function'''
    parser = argparse.ArgumentParser(description="Run REMC simulation.")

    # Required arguments
    parser.add_argument("sequence", type=str, help="The protein sequence in HP format.")
    parser.add_argument("protein_name", type=str, help="The name of the protein.")
    parser.add_argument("optimal_E", type=float, help="The optimal energy level - if unknown a value must still be supplied.")

    # Optional arguments with defaults
    parser.add_argument("--min_t", type=int, default=160, help="Minimum temperature (default: 160).")
    parser.add_argument("--max_t", type=int, default=220, help="Maximum temperature (default: 220).")
    parser.add_argument("--mc_steps", type=int, default=500, help="Monte Carlo steps (default: 500).")
    parser.add_argument("--replica_number", type=int, default=25, help="Number of replicas (default: 25).")
    parser.add_argument("--cutoff_runtime", type=int, default=120, help="Cutoff runtime in seconds (default: 120).")
    parser.add_argument("--nb_of_runs", type=int, default=5, help="Number of runs (default: 5).")

    # Parse the arguments from the command line
    args = parser.parse_args()

    # Call the remc_full_sim function with the parsed arguments
    remc_full_sim(
        sequence=args.sequence,
        protein_name=args.protein_name,
        optimal_E=args.optimal_E,
        min_t=args.min_t,
        max_t=args.max_t,
        mc_steps=args.mc_steps,
        replica_number=args.replica_number,
        cutoff_runtime=args.cutoff_runtime,
        nb_of_runs=args.nb_of_runs
    )


if __name__ == "__main__":
    main()