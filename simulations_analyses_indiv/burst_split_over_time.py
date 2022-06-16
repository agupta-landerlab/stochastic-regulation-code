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

#tracking populations split by whether TF was bursting at t=0

df_use_burstsplit = updated_df

#tracking median across cells that DID vs DID NOT have a burst at t=0 (fig 1D)
t0 = df_use_burstsplit[df_use_burstsplit['sampling_time']==0]

cells_with_burst = list(t0[t0['TF_is_bursting']==True]['Unnamed: 0'])
cells_without_burst_all = list(set(t0['Unnamed: 0']) - set(cells_with_burst))
cells_without_burst = cells_without_burst_all

burst_cells_all_times = []
for j in range(0,40):
    for i in cells_with_burst:
        burst_cells_all_times.extend([int(i+j*20e3)])
        
nonburst_cells_all_times = []
for j in range(0,40):
    for i in cells_without_burst:
        nonburst_cells_all_times.extend([int(i+j*20e3)])
        
print(len(burst_cells_all_times), len(nonburst_cells_all_times))

bursted_cells_info = df_use_burstsplit.loc[burst_cells_all_times,['TF_is_bursting','Target_is_bursting',
                           'total_TF_mRNA','TF_protein_1K','total_Target_mRNA','unspliced_Target','sampling_time']]

nonbursted_cells_info = df_use_burstsplit.loc[nonburst_cells_all_times,['TF_is_bursting','Target_is_bursting',
                           'total_TF_mRNA','TF_protein_1K','total_Target_mRNA','unspliced_Target','sampling_time']]

#TF MRNA
plt.figure(figsize=(18,5))
plt.plot(bursted_cells_info.groupby('sampling_time')['total_TF_mRNA'].mean(),color='darkblue',linewidth=5,label='had burst at t=0')
plt.plot(nonbursted_cells_info.groupby('sampling_time')['total_TF_mRNA'].mean(),color='lightgrey',linewidth=5,label='no burst at t=0')

plt.legend(fontsize=25)
plt.ylabel("# TF mRNA transcripts")
plt.xlabel("time after burst (hrs)")
plt.ylim(0,85)
plt.xticks(ticks=np.arange(0,1300,300),labels=np.arange(0,20,5),fontsize=25)
sns.despine()
plt.show()

#TF PROTEIN
plt.figure(figsize=(18,5))
plt.plot(bursted_cells_info.groupby('sampling_time')['TF_protein_1K'].mean(),color='#6FCAC9',linewidth=5,label='had burst at t=0')
plt.plot(nonbursted_cells_info.groupby('sampling_time')['TF_protein_1K'].mean(),color='lightgrey',linewidth=5,label='no burst at t=0')

plt.legend(fontsize=25)
plt.ylabel("# TF protein molecules (1k)")
plt.ylim(30,115)
plt.xticks(ticks=np.arange(0,1300,300),labels=np.arange(0,20,5),fontsize=25)
plt.xlabel("time after burst (hrs)")
sns.despine()
plt.show()

#TARGET MRNA
plt.figure(figsize=(18,5))
plt.plot(bursted_cells_info.groupby('sampling_time')['total_Target_mRNA'].mean(),color='red',linewidth=5,label='had burst at t=0')
plt.plot(nonbursted_cells_info.groupby('sampling_time')['total_Target_mRNA'].mean(),color='lightgrey',linewidth=5,label='no burst at t=0')

plt.plot(bursted_cells_info.groupby('sampling_time')['unspliced_Target'].mean(),color='red',linewidth=5,label='had burst at t=0')
plt.plot(nonbursted_cells_info.groupby('sampling_time')['unspliced_Target'].mean(),color='lightgrey',linewidth=5,label='no burst at t=0')


plt.legend(fontsize=25)
plt.ylabel("# Target mRNA transcripts")
plt.ylim(0,85)
plt.xticks(ticks=np.arange(0,1300,300),labels=np.arange(0,20,5),fontsize=25)
plt.xlabel("time after burst (hrs)",fontsize=22)
sns.despine()
plt.show()