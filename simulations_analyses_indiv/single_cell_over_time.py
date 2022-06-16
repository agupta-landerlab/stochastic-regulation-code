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

#tracking the same cell over time

random_cells_check = [176,170,15,47]
for i in random_cells_check:
    print(i)
    rows_keep = list(np.arange(i,len(updated_df),20e3))
    single_cell_over_time = updated_df.loc[rows_keep,['TF_is_bursting','Target_is_bursting','total_TF_mRNA','TF_protein_1K','total_Target_mRNA','sampling_time']]

    #TF mRNA
    f2b_TFRNA = plt.figure(figsize=(10,2))
    plt.plot(single_cell_over_time['sampling_time'], single_cell_over_time['total_TF_mRNA'],color='darkblue',linewidth=3)
    plt.ylabel("TF mRNA\nabundance")
    plt.xticks(np.arange(0,1500,300),np.arange(0,25,5),fontsize=18)
    plt.xlabel("time (hour)")
    sns.despine()
    plt.show()

    #TF protein
    f2b_TFprot = plt.figure(figsize=(10,2))
    plt.plot(single_cell_over_time['sampling_time'], single_cell_over_time['TF_protein_1K'],
             color='teal',linewidth=4,alpha=0.5)
    plt.ylabel("#TF proteins\n(1k molecules)")
    plt.xticks(np.arange(0,1500,300),np.arange(0,25,5),fontsize=18)
    plt.xlabel("time (hour)")
    sns.despine()
    plt.show()

    #Target mRNA
    f2b_TargetRNA = plt.figure(figsize=(10,2))
    plt.plot(single_cell_over_time['sampling_time'], single_cell_over_time['total_Target_mRNA'],color='red',linewidth=3)
    plt.ylabel("Target\n mRNA abundance")
    plt.xticks(np.arange(0,1500,300),np.arange(0,25,5),fontsize=18)
    plt.xlabel("time (hour)")
    sns.despine()
    plt.show()