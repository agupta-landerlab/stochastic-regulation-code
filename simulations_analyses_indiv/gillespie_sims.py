from gillespie_vCaltech_dependencies import *
import argparse
import numpy as np

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--n_cells', type=int)
parser.add_argument('--n_mins', type=int)
parser.add_argument('--batch', type=int)
args = parser.parse_args()
print('\n\n****')
print(args)
print('****\n\n')

# Column 0 is change in TF mRNA (m), column 1 is change in TF protein (p), column 2 is change in Target mRNA (m2)
# Each row represents a reaction
simple_update = np.array(
    [
        [1, 0, 0, 0, 0], # turn TF txn ON, if in OFF state
        [-1, 0, 0, 0, 0], # turn TF txn OFF, if in ON state
        [0, 1, 0, 0, 0],  # Make TF mRNA transcript
        [0, -1, 0, 0, 0],  # Degrade TF mRNA
        [0, 0, 1, 0, 0],  # Make TF protein
        [0, 0, -1, 0, 0],  # Degrade TF protein
        [0, 0, 0, 1, 0], # Make Target mRNA transcript
        [0, 0, 0, -1, 0], # Degrade Target mRNA
        [0, 0, 0, 0, 1], # turn Target txn ON, if in OFF state
        [0, 0, 0, 0, -1] # turn Target txn OFF, if in ON state
    ], dtype=np.int)

# Specify parameters for calculation (use median values from literature)

## to incorporate burst switching
# TF
kon = 0.27*60
koff = 8.4*60
burst_size = 32
tsn_rate = 0.059*60
# Target
kon2 = 0.25*60
koff2 = 7.7*60
burst_size2 = 40

# TF RNA, TF protein, Target RNA production rates (txn = k_off*burst_size; tsn rate = protein production rate)
beta_m = koff * burst_size
beta_p = tsn_rate
beta_m2 = koff2 * burst_size2

# TF RNA, TF protein, Target RNA decay rates (using median half-lives)
gamma_m = np.log(2)/(2.5/60)
gamma_p = np.log(2)/(24/60)
gamma_m2 = np.log(2)/(3.7/60)

sim_args = (kon, koff, kon2, koff2, beta_m, beta_p, beta_m2, gamma_m, gamma_p, gamma_m2) # txn rate, tsn rate, decay rate
time_points = np.linspace(0, args.n_mins, int(args.n_mins/10)) ## UPDATE AS NEEDED
population_0 = np.array([0, 0, 0, 0, 0], dtype=int) #TF bursting state, TF RNA, TF protein, Target RNA
size = args.n_cells ## UPDATE

# Seed random number generator for reproducibility
np.random.seed(42)

# Initialize output array
samples = np.empty((size, len(time_points), 5), dtype=int)

# Run the calculations
for i in tqdm.tqdm(range(size)):
    samples[i, :, :] = gillespie_ssa(simple_propensity, simple_update, population_0, time_points, args=sim_args)
    
import pandas as pd

all_samples = pd.DataFrame()
TF_RNA, TF_P, TARGET_RNA, TIME = [],[],[],[]
for tp in range(len(time_points)):
    this_tp_samples = samples[:, tp, :]
    TF_RNA.extend(this_tp_samples[:,1])
    TF_P.extend(this_tp_samples[:,2])
    TARGET_RNA.extend(this_tp_samples[:,3])
    TIME.extend([tp]*len(this_tp_samples))
    
all_samples['TF_RNA'] = TF_RNA
all_samples['TF_PROTEIN'] = TF_P
all_samples['Target_RNA'] = TARGET_RNA
all_samples['tp'] = TIME

#code to reshape 3D to 2D and output + input: geeksforgeeks.org/how-to-load-and-save-3d-numpy-array-to-file-using-savetxt-and-loadtxt-functions/
samples_reshaped = samples.reshape(samples.shape[0], -1)

#outputting the np array
np.savetxt('/data/srlab/agupta/lander_lab/gillespie_outputs/NP_samples_20h_'+str(args.batch)+'.csv',
                samples_reshaped, delimiter=",")

#outputting the DF
all_samples.to_csv("/data/srlab/agupta/lander_lab/gillespie_outputs/DF_samples_20h_"+str(args.batch)+".csv")