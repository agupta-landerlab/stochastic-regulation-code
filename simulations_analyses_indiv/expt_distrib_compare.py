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

# comparing simulated TF RNA, TF protein, and Target RNA median abundances to thousands of literature-based genes

df_use_S1 = updated_df

sim_data_TF_t0 = df_use_S1[df_use_S1['sampling_time']==0]['total_TF_mRNA'][::200]
sim_data_Target_t0 = df_use_S1[df_use_S1['sampling_time']==0]['total_Target_mRNA']
sim_data_TF_PROT_t0 = [i for i in df_use_S1[df_use_S1['sampling_time']==0]['TF_protein_1K']]

#data from Schwanhausser et al
nature_2011 = pd.read_excel('parameters_extraction/41586_2011_BFnature10098_MOESM304_ESM.xls')
real_PROT = [i for i in nature_2011['Protein copy number average [molecules/cell]']]
real_RNA = nature_2011['mRNA copy number average [molecules/cell]']

print(np.median(sim_data_TF_PROT_t0),np.median([i/1000 for i in real_PROT]))

plt.figure(figsize=(5,5))
plt.scatter(np.arange(0,4309),sorted(np.log2(real_RNA.dropna())),color='black',alpha=0.2,s=2)
plt.scatter(2800,np.log2(np.median(sim_data_TF_t0)),color='darkblue',s=70)
plt.scatter(3350,np.log2(np.median(sim_data_Target_t0)),color='red',s=70)
plt.xlabel("genes sorted by mean\nmRNA abundance")
plt.ylabel("log2(mRNA transcripts)")
sns.despine()
plt.show()

plt.figure(figsize=(5,5))
plt.scatter(np.arange(0,5000),sorted([np.log10(i/1000) for i in real_PROT]),color='black',alpha=0.2,s=2)
plt.scatter(2930,np.median(np.log10(sim_data_TF_PROT_t0)),color='lightblue',s=70)
plt.xlabel("genes sorted by mean\nprotein abundance")
plt.ylabel("log10(protein abundance/1k)")
sns.despine()
plt.show()

#what percentile of genes are represented by the simulated genes
#TF MRNA, TF PROTEIN, TARGET MRNA
2800/4309, 2930/5000, 3350/4309


# comparing TF RNA and Target RNA steady-state distributions in simulations to experimentally-derived values for an example TF and non-TF gene

df_use_QQ = updated_df

nanog_smFISH = pd.read_csv("parameters_extraction/Nanog_smFISH.csv")
#col names: NM_LIF	NM_2i	B6_LIF	B6_2i
nanog_real_counts = nanog_smFISH['B6_LIF'].dropna()

plt.figure(figsize=(6,3))
plt.plot([0,180],[0,180],color='red',linewidth=2)
plt.scatter(sorted(nanog_real_counts),sorted(sim_data_TF_t0),color='black')
plt.xlabel("real counts")
plt.ylabel("simulated counts")
sns.despine()
# plt.title("TF\n(NANOG)")
plt.show()

TNFR1_df = pd.DataFrame(pd.read_csv("parameters_extraction/TNFR1_mRNA_prot.csv"))
sim_data_Target_t0_subset = sim_data_Target_t0[::540]
plt.figure(figsize=(6,3))
plt.plot([0,200],[0,200],color='red',linewidth=2)
plt.scatter(sorted(TNFR1_df['RNA']),sorted(sim_data_Target_t0_subset),color='black')
plt.xlabel("real counts")
plt.ylabel("simulated counts")
sns.despine()
plt.xlim(0,200)
# plt.title("non-TF\n(TNFR1)")