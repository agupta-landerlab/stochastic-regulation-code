import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import random
import json
from stochastic_bursting_simulator import simulate
from multiprocessing import Pool

import multiprocessing
from time import sleep

#default parameter values are median from literature
params = {'k_on_TF': 0.27,
          'k_off_TF': 8.4,
          'burst_size_TF': 32,
          'k_on_Target': 0.25,
          'k_off_Target': 7.7,
          'burst_size_Target': 40,
          'splicing_half_life_minutes': 7,
          'mrna_half_life_TF': 2.5,
          'mrna_half_life_Target': 3.7,
          'protein_half_life': 28, 
          'protein_production_rate': 0.059,
          'labeling_efficiency': 0.8,
          'pulse_time': 60,
          'num_cells': 20_000,
          'dynamics': 'MM',
          'capture_efficiency': 1}

def caller(args_list):
    args, filename = args_list
    all_samples = simulate(**args)
    all_samples.to_csv(filename+".csv")

if __name__ == '__main__':
    #optional: run sim many times
    num_runs = 25
    args_list = []
    this_param_quartile_args = params.copy()
    for sim_run in range(num_runs):
	#update based on param settings
        filename = "sim_outputs/samples_over_time_"+str(sim_run)
        args_list.append((this_param_quartile_args, filename))
                                 
# loop over list of parameters
    with multiprocessing.Pool(processes=32) as pool:
        pool.map(caller, args_list)
