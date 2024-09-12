'''Main script to run to execute REMC simulation
'''

# Import of necessary modules

import numpy as np
import math as m
import os
import sys

import copy
import random
import time
import multiprocessing as mp
from tqdm import tqdm
import matplotlib.pyplot as plt

from components import *

# Main program function

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



if __name__ == "__main__":
    print('gsjdfl')