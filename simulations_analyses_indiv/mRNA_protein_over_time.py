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

#JOINT DISTRIBUTIONS & CORRELATIONS

pulse_length = 60

def joint_distribs(a):
    end, start = list(a[a['sampling_time']==pulse_length]['mRNA_ever_produced_TF']), list(a[a['sampling_time']==0]['mRNA_ever_produced_TF'])
    new_during_pulse = [end[i]-start[i] for i in range(len(end))]

    tot_mRNA_prot_corrs, new_mRNA_prot_corrs = [], []
    mRNA_autocorrs, prot_autocorrs = [], []

    for tp in np.arange(0,1440,30):
        #TF mRNA:protein correlations (total, ∆)
        tot_TF_mRNA = a[a['sampling_time']==0]['total_TF_mRNA']
        tot_TF_prot = a[a['sampling_time']==tp]['TF_protein_1K']
        
        tot_mRNA_prot_corr = spearmanr(tot_TF_mRNA, tot_TF_prot)
        new_mRNA_prot_corr = spearmanr(new_during_pulse, tot_TF_prot)

        #autocorrelations (TF mRNA:TF mRNA over time & TF protein:TF protein over time)
        mRNA_autocorr = spearmanr(a[a['sampling_time']==tp]['total_TF_mRNA'],tot_TF_mRNA)
        prot_autocorr = spearmanr(a[a['sampling_time']==tp]['TF_protein_1K'],a[a['sampling_time']==0]['TF_protein_1K'])

        tot_mRNA_prot_corrs.append(tot_mRNA_prot_corr[0])
        new_mRNA_prot_corrs.append(new_mRNA_prot_corr[0])
        mRNA_autocorrs.append(mRNA_autocorr[0])
        prot_autocorrs.append(prot_autocorr[0])   
        
        if tp==420:
            plt.figure(figsize=(8,8))
            corr = spearmanr(tot_TF_mRNA, tot_TF_prot)[0]
            sns.jointplot(tot_TF_mRNA, tot_TF_prot, color='#6F63AC',alpha=0.5)
            sns.despine()
            plt.show()
    
    return tot_mRNA_prot_corrs, new_mRNA_prot_corrs, mRNA_autocorrs, prot_autocorrs

## tracking mRNA and protein correlations over time (one and many sim runs)

#ONE RUN
tot_mRNA_prot, new_mRNA_prot, mRNA_auto, prot_auto = joint_distribs(updated_df)

plt.figure(figsize=(12,8))
plt.plot(mRNA_auto, label = 'TF mRNA autocorrelation',color='darkblue',linewidth=7)
plt.plot(prot_auto, label = 'TF protein autocorrelation',color='#6FCAC9',linewidth=7)
plt.plot(tot_mRNA_prot, label = 'total TF mRNA:total TF protein',color='mediumpurple',linewidth=7)
plt.xticks(np.arange(0,48,8),np.arange(0,24,4))
sns.despine()
plt.legend(fontsize=25)
plt.ylabel("Spearman ρ",fontsize=28)
plt.xlabel("Hours after first entity is measured",fontsize=25)
plt.show()

#MANY (25) RUNS
# for each corr type, create lineplot
corr_names = ['TF_mRNA:protein','RNA:RNA','protein:protein','TF_total_t0:Target_total_tN','TF_total_t0:Target_nascent_tN','Target_total_t0:TF_total_tN','Target_total_t0:TF_nascent_tN']
    
all_runs_df = pd.DataFrame()
run_nums=[]
for sim_run in range(0,25):
    filename = 'simulated_data_for_figures/fig2_and_3_medians/medians_many_runs/spearman_over_time_run'+str(sim_run)+'.csv'
    this_run_corrs_df = pd.DataFrame(pd.read_csv(filename))
    
    all_runs_df = pd.concat([all_runs_df, this_run_corrs_df])
    run_nums.extend([sim_run]*len(this_run_corrs_df))
    
all_runs_df['run_num'] = run_nums
all_runs_df['timepoint']=[int(i/30) for i in all_runs_df['timepoint']]
    
f = plt.figure(figsize=(12,8))
sns.lineplot(data=all_runs_df,x='timepoint',y='TF_mRNA:protein',color='mediumpurple', ci=95,lw=2)#, err_style='bars')
sns.lineplot(data=all_runs_df,x='timepoint',y='RNA:RNA',color='darkblue', ci=95,lw=2)#, err_style='bars')
sns.lineplot(data=all_runs_df,x='timepoint',y='protein:protein',color='#6FCAC9', ci=95,lw=2)#,err_style='bars')

sns.stripplot(data=all_runs_df,x='timepoint',y='TF_mRNA:protein',color='mediumpurple',s=4)
sns.stripplot(data=all_runs_df,x='timepoint',y='RNA:RNA',color='darkblue',s=4)
sns.stripplot(data=all_runs_df,x='timepoint',y='protein:protein',color='#6fcac9',s=4)
plt.xticks(np.arange(0,48,8),np.arange(0,24,4))
sns.despine()
plt.legend(fontsize=25)
plt.ylabel("Spearman ρ",fontsize=28)
plt.xlabel("Hours after first entity is measured",fontsize=25)