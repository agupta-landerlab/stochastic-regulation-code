from sim_dependencies.py import *

#read in simulated data

sim_file_use = 'simulated_data_for_figures_fig2_and_3_medians_20k_30min_poisson_decay_SAMPLES_OVER_TIME.csv'
num_sampling_times = 48
num_cells = 20e3

updated_df = pd.read_csv(sim_file_use)
updated_df['unspliced_Target'] = updated_df['unspliced_unlabeled_Target']+updated_df['unspliced_labeled_Target']
updated_df['unspliced_TF'] = updated_df['unspliced_unlabeled_TF']+updated_df['unspliced_labeled_TF']

updated_df['labeled_Target'] = updated_df['unspliced_labeled_Target']+updated_df['spliced_labeled_Target']
updated_df['labeled_TF'] = updated_df['unspliced_labeled_TF']+updated_df['spliced_labeled_TF']

updated_df['unlabeled_Target'] = updated_df['unspliced_unlabeled_Target']+updated_df['spliced_unlabeled_Target']
updated_df['unlabeled_TF'] = updated_df['unspliced_unlabeled_TF']+updated_df['spliced_unlabeled_TF']

#effect of varying gene-specific parameters on steady-state distributions for the 3 key entities
#combine DFs with appropriate values, from each quartile

for param in ["protein_half_life", "mrna_half_life_TF", "protein_production_rate",
              "k_on_TF", "k_off_TF", "burst_size_TF", "burst_size_Target", "k_on_Target",
              "k_off_Target", "mrna_half_life_Target"]:
    full_param_df = pd.DataFrame()

    for quartile in ['Q1','med','Q3']:
        steady_state_distribs_file = "STEADY_STATE_DISTRIB_5k_cells_"+param+"_"+quartile+".csv"
        this_quartile_df = pd.read_csv("simulated_data_for_figures/figS1/"+steady_state_distribs_file)
        S2_num_cells = len(this_quartile_df)
        this_quartile_df = pd.DataFrame(pd.concat([this_quartile_df['total_TF_mRNA'],
                                      this_quartile_df['TF_protein_1K'],
                                      this_quartile_df['total_Target_mRNA']]))
        full_param_df = pd.concat([full_param_df, this_quartile_df])
    
    full_param_df['entity'] = np.reshape((['total TF mRNA']*S2_num_cells + ['TF protein 1K']*S2_num_cells + ['total Target mRNA']*S2_num_cells)*3, -1, 1)
    full_param_df['quartile'] = np.reshape(['Q1']*S2_num_cells*3+['median']*S2_num_cells*3+['Q3']*S2_num_cells*3,-1,1)
    
    plt.figure(figsize=(6,6))
    plt.title("Varying "+ " ".join(param.split("_")), fontsize=25)
    plt.xlabel("Entity", fontsize=20)
    sns.boxplot(data=full_param_df, x='entity', y=0, hue='quartile', palette="Greys")
    plt.xticks(rotation=15)
    plt.xlabel("")
    plt.ylabel("Abundances")
    sns.despine()
    plt.legend(loc='upper left')
    plt.show()
    
#steady-state distributions at one timepoint (using median parameter values)

plt.figure(figsize=(8,10))
for t in [0]:
    med_df_t0 = updated_df[updated_df['sampling_time']==t]
    sns.distplot(med_df_t0['total_Target_mRNA'],color='Red',
                 bins=100,hist=False,kde=True,
                 kde_kws={"color": "Red", "lw": 7})
    sns.distplot(med_df_t0['total_TF_mRNA'],color='DarkBlue',
                 bins=100,hist=False,kde=True,
                 kde_kws={"color": "DarkBlue", "lw": 7})
    sns.distplot(med_df_t0['TF_protein_1K'], color='#6FCAC9',
                 bins=100,hist=False,kde=True,
                 kde_kws={"lw": 7})
    sns.despine()
    plt.xlabel("abundance")
    plt.ylabel("density\n(fraction of cells)")
    plt.xlim(0,270)
    plt.ylim(0,0.025)
plt.show()