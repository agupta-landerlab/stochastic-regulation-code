from sim_dependencies.py import *

# read in simulated data from many runs with each parameter value setting
# note k off file Q1 and Q3 have been swapped to reflect burst durations instead!

# fold-change (interdecile ratio D10:D1)
for param in ["protein_half_life", "mrna_half_life_TF", "protein_production_rate",
              "k_on_TF", "k_off_TF", "burst_size_TF", "burst_size_Target", "k_on_Target",
              "k_off_Target", "mrna_half_life_Target"]:

    all_quartiles_df = pd.DataFrame()
    for quartile in ['Q1','med','Q3']:
        
        param_df = pd.read_csv("simulated_data_for_figures/fig4_burstparams/10k_cells_"+param+"_"+quartile+"_ES_over_time.csv")
        num_timepoints = len(set(param_df['timepoint']))
        sim_runs = []
        for i in range(10):
            sim_runs.extend([i]*num_timepoints)
        param_df['sim_run'] = sim_runs
        param_df['quartile'] = [quartile]*len(param_df)
         t
        all_quartiles_df = pd.concat([all_quartiles_df, param_df])

    plt.figure(figsize=(6,6))
    sns.boxplot(data=all_quartiles_df,x='quartile',
                y='Estimated Effect Size', 
                palette=['#7FC993','#6F63AC','#3C81A3'],notch=True)
    plt.ylabel("Fold Change")
    plt.xlabel("")
    sns.despine()
    plt.title("Varying "+param)
    plt.show()
    
# Spearman rho over time

for param in ["mrna_half_life_TF", "protein_half_life", 
              "k_on_TF", "k_off_TF", "k_on_Target", "mrna_half_life_Target","k_off_Target"]:

    all_quartiles_df = pd.DataFrame()
    max_vals, max_ids = [],[]
    for quartile in ['Q1','med','Q3']:
        param_df = pd.read_csv("simulated_data_for_figures/fig4_burstparams/20k_cells_"+param+"_"+quartile+"_Spearman_over_time.csv")
        num_timepoints = len(set(param_df['timepoint']))
        sim_runs = []
        for i in range(25):
            sim_runs.extend([i]*num_timepoints)
        param_df['sim_run'] = [int(i) for i in sim_runs]
        param_df['quartile'] = [quartile]*len(param_df)
        
        cols_keep = ['TF_total_t0:Target_total_tN','TF_total_t0:Target_nascent_tN','timepoint','sim_run','quartile']
        param_df = param_df[cols_keep]
        
        max_vals.append(np.max(param_df['TF_total_t0:Target_nascent_tN'])) #param_df.groupby('timepoint').mean()['TF_total_t0:Target_nascent_tN']))#
        max_ids.append(np.argmax(param_df['TF_total_t0:Target_nascent_tN'])) #param_df.groupby('timepoint').mean()['TF_total_t0:Target_nascent_tN']))#
        
        all_quartiles_df = pd.concat([all_quartiles_df, param_df])
    
    f4line_params = plt.figure(figsize=(10,6))
    ax1 = sns.lineplot(data=all_quartiles_df,
                 x='timepoint',
                 y='TF_total_t0:Target_nascent_tN', hue='quartile',
                 palette=['#7fc993','#6f63ac','#3c81a3'],linewidth=3,err_style='bars')

    legend1=ax1.legend()
    legend1.remove()
    plt.xticks(np.arange(0,1440,240),np.arange(0,24,4))
    plt.ylabel("Spearman ρ")
    plt.xlabel("")
    sns.despine()
    plt.ylim(0,0.3)
    plt.title("Varying "+param)
    plt.show()
    
    f4jitter_params = plt.figure(figsize=(10,6))
    ax2 = sns.stripplot(data=all_quartiles_df,x='timepoint',y='TF_total_t0:Target_nascent_tN',
                  hue='quartile',palette=['#3c81a3','#6f63ac','#7fc993'],alpha=0.5)
    legend2=ax2.legend()
    legend2.remove()
    plt.ylabel("Spearman ρ")
    plt.xlabel("")
    sns.despine()
    plt.title("Varying "+param)
    plt.show()