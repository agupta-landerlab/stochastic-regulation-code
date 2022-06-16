from sim_dependencies.py import *

#run simulator 2X, with NO link function, and with different basal TF and Target abundances
#(by varying k_on for both TF and Target from Q1 to Q3)
folder = 'simulated_data_for_figures/fig4C_states/many_runs/'

all_state_based_df = pd.DataFrame()

i=0
for Q1_fn, Q3_fn in [("k_on_Q1_STATE_and_REG_samplesOverTime.csv","k_on_Q3_STATE_and_REG_samplesOverTime.csv"),
                     ('k_on_Q1_SAMPLES_OVER_TIME_MANY_SAMPLES.csv','k_on_Q3_SAMPLES_OVER_TIME_MANY_SAMPLES.csv')]:

    ####### READ IN AND UPDATE HIGH + LOW DFS ########
    print('reading in DFs for this version of state-based')
    states_low_df = pd.DataFrame(pd.read_csv(folder+Q1_fn))
    states_high_df = pd.DataFrame(pd.read_csv(folder+Q3_fn))

    states_low_df['state'] = ['low']*len(states_low_df)
    states_low_df['unspliced_Target'] = states_low_df['unspliced_unlabeled_Target']+states_low_df['unspliced_labeled_Target']
    states_low_df['unspliced_TF'] = states_low_df['unspliced_unlabeled_TF']+states_low_df['unspliced_labeled_TF']
    states_low_df['unlabeled_Target'] = states_low_df['unspliced_unlabeled_Target']+states_low_df['spliced_unlabeled_Target']
    states_low_df['labeled_TF'] = states_low_df['unspliced_labeled_TF']+states_low_df['spliced_labeled_TF']
    states_low_df['unlabeled_TF'] = states_low_df['unspliced_unlabeled_TF']+states_low_df['spliced_unlabeled_TF']
    states_low_df['labeled_Target'] = states_low_df['unspliced_labeled_Target']+states_low_df['spliced_labeled_Target']

    states_high_df['state'] = ['high']*len(states_high_df)
    states_high_df['unspliced_Target'] = states_high_df['unspliced_unlabeled_Target']+states_high_df['unspliced_labeled_Target']
    states_high_df['unspliced_TF'] = states_high_df['unspliced_unlabeled_TF']+states_high_df['unspliced_labeled_TF']
    states_high_df['unlabeled_Target'] = states_high_df['unspliced_unlabeled_Target']+states_high_df['spliced_unlabeled_Target']
    states_high_df['labeled_TF'] = states_high_df['unspliced_labeled_TF']+states_high_df['spliced_labeled_TF']
    states_high_df['unlabeled_TF'] = states_high_df['unspliced_unlabeled_TF']+states_high_df['spliced_unlabeled_TF']
    states_high_df['labeled_Target'] = states_high_df['unspliced_labeled_Target']+states_high_df['spliced_labeled_Target']

    ### POISSON DOWNSAMPLE COUNTS ###
    print('downsampling counts')
    DS_FRAC = 0.1
    for col in ['total_TF_mRNA','total_Target_mRNA','labeled_TF','unlabeled_Target',
               'unlabeled_TF','labeled_Target','unspliced_TF','unspliced_Target']:
        states_high_df[col] = poisson_capture(DS_FRAC,states_high_df[col])
        states_low_df[col] = poisson_capture(DS_FRAC,states_low_df[col])

    #############################################
    all_runs_state_corrs_df = pd.DataFrame()
    run_nums=[]
    print('running through sampling runs')
    for sampling_run in range(0,25):
        num_cells_per = 10000
        print('subsetting from high and low DFs, to create 2 "states"')
        low_rand_cells_i_t0 = np.arange(sampling_run*num_cells_per,
                                        sampling_run*num_cells_per+num_cells_per)
        high_rand_cells_i_t0 = np.arange(sampling_run*num_cells_per,
                                         sampling_run*num_cells_per+num_cells_per)

        #track cell indices at all timepoints, ST tracking the same cells as t0
        low_cells_all_tps, high_cells_all_tps = [], []
        for TP in range(0,48): #tps (for 1440 min, w 30 min intervals)
            total_cells_per_tp = 25e4
            low_cells_all_tps.extend(list(low_rand_cells_i_t0+TP*total_cells_per_tp))
            high_cells_all_tps.extend(list(high_rand_cells_i_t0+TP*total_cells_per_tp))
        
        states_low_df_sampled = states_low_df.iloc[low_cells_all_tps,:]
        states_high_df_sampled = states_high_df.iloc[high_cells_all_tps,:]

        print('combining into a joint df')
        #combine cells from high and low states at each timepoint
        joint_states_df = pd.concat([states_low_df_sampled,states_high_df_sampled])

        # FIGURE 5C
        tp_check = sorted(list(set(joint_states_df['sampling_time'])))
        joint_states_df_t0 = joint_states_df[joint_states_df['sampling_time']==0]
        print('calculating corrs over time')
        total_corrs=[]
        new_corrs=[]
        REVERSE_new_corrs=[]
        LU_corrs,LT_corrs,UL_corrs,TL_corrs,TT_corrs=[],[],[],[],[]
        for tp in tp_check:
            tp_data = joint_states_df[joint_states_df['sampling_time']==tp]
            Target_tp_data = tp_data['total_Target_mRNA']
            TF_tp_data = tp_data['total_TF_mRNA']
            new_Target_tp_data = tp_data['unspliced_Target']
            new_TF_tp_data = tp_data['unspliced_TF']
            labeled_TF_data = tp_data['labeled_TF']
            unlabeled_Target_data = tp_data['unlabeled_Target']
            unlabeled_TF_data = tp_data['unlabeled_TF']
            labeled_Target_data = tp_data['labeled_Target']

            total_corrs.append(spearmanr(joint_states_df_t0['total_TF_mRNA'],Target_tp_data)[0])
            new_corrs.append(spearmanr(joint_states_df_t0['total_TF_mRNA'],new_Target_tp_data)[0])
            REVERSE_new_corrs.append(spearmanr(joint_states_df_t0['total_Target_mRNA'],new_TF_tp_data)[0])   

            LU_corrs.append(spearmanr(labeled_TF_data,unlabeled_Target_data)[0])
            LT_corrs.append(spearmanr(labeled_TF_data,Target_tp_data)[0])
            UL_corrs.append(spearmanr(unlabeled_TF_data,labeled_Target_data)[0])
            TL_corrs.append(spearmanr(TF_tp_data,labeled_Target_data)[0])
            TT_corrs.append(spearmanr(TF_tp_data,Target_tp_data)[0])

        sampling_times = tp_check
        if i==1:
            this_run_corrs_df = pd.DataFrame(list(zip(new_corrs, REVERSE_new_corrs, sampling_times)),
                                             columns=['state-new-corrs','state-reverse-corrs','sampling_time'])
        else:
            this_run_corrs_df = pd.DataFrame(list(zip(new_corrs, REVERSE_new_corrs, sampling_times)),
                                             columns=['state-reg-new-corrs','state-reg-reverse-corrs','sampling_time'])

        run_nums.extend([sampling_run]*len(this_run_corrs_df))

        all_runs_state_corrs_df = pd.concat([all_runs_state_corrs_df, this_run_corrs_df])

        all_runs_state_corrs_df['sim_run'] = run_nums

    all_runs_state_corrs_df['sampling_time']=[int(i/30) for i in all_runs_state_corrs_df['sampling_time']]
    all_state_based_df = pd.concat([all_state_based_df, all_runs_state_corrs_df],axis=1)
    i+=1

    
# all_state_based_df.to_csv("simulated_data_for_figures/fig4C_states/many_runs/states_combined_downsampled_CORRS_df.csv")
all_state_based_df = pd.DataFrame(pd.read_csv("simulated_data_for_figures/fig4C_states/many_runs/states_combined_downsampled_CORRS_df.csv"))

#### REGULATORY IMPORT ####
#get forward and reverse corrs from DOWNSAMPLED counts!! 10% for now
# corr_names = ['TF_mRNA:protein','RNA:RNA','protein:protein','TF_total_t0:Target_total_tN',
#               'TF_total_t0:Target_nascent_tN','Target_total_t0:TF_total_tN','Target_total_t0:TF_nascent_tN']
    
# all_runs_DOWNSAMPLED_df = pd.DataFrame()
# run_nums=[]
# for sim_run in range(0,25):
#     filename = 'simulated_data_for_figures/fig4C_downsampling/old/rhos_capture_eff=0.1_run'+str(sim_run)+'.csv'
#     this_run_corrs_df = pd.DataFrame(pd.read_csv(filename))
    
#     all_runs_DOWNSAMPLED_df = pd.concat([all_runs_DOWNSAMPLED_df, this_run_corrs_df])
#     run_nums.extend([sim_run]*len(this_run_corrs_df))
    
# all_runs_DOWNSAMPLED_df['run_num'] = run_nums
# all_runs_DOWNSAMPLED_df['timepoint']=[int(i/30) for i in all_runs_DOWNSAMPLED_df['timepoint']]

fig4_states = plt.figure(figsize=(12,8))

####### REGULATORY #######
#FORWARD, downsampled
sns.lineplot(data=all_runs_DOWNSAMPLED_df,x='timepoint',y='TF_total_t0:Target_nascent_tN',color='#6f63ac',
             ci=95,lw=2, err_style='bars',label='forward regulatory')
sns.stripplot(data=all_runs_DOWNSAMPLED_df,x='timepoint',y='TF_total_t0:Target_nascent_tN',color='#6f63ac',s=4)

#REVERSE
sns.lineplot(data=all_runs_DOWNSAMPLED_df,x='timepoint',y='Target_total_t0:TF_nascent_tN',color='grey',
             ci=95,lw=2, err_style='bars',label='reverse regulatory')
sns.stripplot(data=all_runs_DOWNSAMPLED_df,x='timepoint',y='Target_total_t0:TF_nascent_tN',color='grey',s=4)

####### STATE-BASED #######
sns.lineplot(data=all_state_based_df,x='sampling_time',y='state-new-corrs',color='#6f63ac',
             ci=95,lw=2,alpha=0.5,dashes=True, err_style='bars',label='forward state only')
sns.lineplot(data=all_state_based_df,x='sampling_time',y='state-reverse-corrs',color='grey',
             ci=95,lw=2,alpha=0.5,dashes=True, err_style='bars',label='reverse state only')
sns.stripplot(data=all_state_based_df,x='sampling_time',y='state-new-corrs',color='#6f63ac',s=4,alpha=0.2)
sns.stripplot(data=all_state_based_df,x='sampling_time',y='state-reverse-corrs',color='grey',s=4,alpha=0.2)

####### STATE-BASED WITH REGULATION ALSO ########
sns.lineplot(data=all_state_based_df,x='sampling_time',y='state-reg-new-corrs',color='red',
             ci=95,lw=2,alpha=.8,dashes=True, err_style='bars',label='forward state + reg')
sns.lineplot(data=all_state_based_df,x='sampling_time',y='state-reg-reverse-corrs',color='salmon',
             ci=95,lw=2,alpha=0.8,dashes=True, err_style='bars',label='reverse state + reg')
sns.stripplot(data=all_state_based_df,x='sampling_time',y='state-reg-new-corrs',color='red',s=4,alpha=0.8)
sns.stripplot(data=all_state_based_df,x='sampling_time',y='state-reg-reverse-corrs',color='salmon',s=4,alpha=0.8)

sns.despine()
plt.legend(fontsize=12,loc='best')
plt.xticks(np.arange(0,48,8),np.arange(0,24,4),fontsize=18)
plt.xlabel("Time (hours)")
plt.ylabel("Spearman ρ")
plt.show()


## cell state boxplots
distrib_boxplot_df = pd.DataFrame()
distrib_boxplot_df['total mRNA'] = list(joint_states_df_t0['total_TF_mRNA'])+list(joint_states_df_t0['total_Target_mRNA'])
distrib_boxplot_df['gene'] = ['TF']*len(joint_states_df_t0) + ['Target']*len(joint_states_df_t0)
distrib_boxplot_df['state'] = ['low']*1000 + ['high']*1000 + ['low']*1000 + ['high']*1000
distrib_boxplot_df

plt.figure(figsize=(5,7))
my_pal = {"TF": "darkblue", "Target": "red"}
sns.violinplot(data=distrib_boxplot_df, x='state', hue='gene', y='total mRNA',palette=my_pal)
plt.legend(loc='upper left')
sns.despine()
plt.show()

from sklearn.decomposition import PCA
pca = PCA(n_components=2)
joint_states_df_t0['total_RNA'] = joint_states_df_t0['total_TF_mRNA']+joint_states_df_t0['total_Target_mRNA']
pca.fit(joint_states_df_t0[['total_TF_mRNA','total_Target_mRNA']])
print(pca.explained_variance_ratio_)
print(pca.singular_values_)
to_transform = joint_states_df_t0.copy()
low_D = pca.fit_transform(to_transform[['total_TF_mRNA','total_Target_mRNA']])

low_D_df = pd.DataFrame(np.array(list([low_D[:,0],low_D[:,1],joint_states_df_t0['total_RNA'],
                                      joint_states_df_t0['total_TF_mRNA'],
                                      joint_states_df_t0['total_Target_mRNA'],
                                      joint_states_df_t0['state']])).T,
                       columns=['PC1','PC2','total_RNA','total_TF_mRNA','total_Target_mRNA','state'])

#for coloring by state
state_colors = [0 if state=='low' else 1 for state in low_D_df['state']]
low_D_df['state color'] = state_colors

low_D_df.plot.scatter(x='PC1',y='PC2', c='total_RNA', cmap='Greys', alpha=0.2, s=20)#'total_RNA'#'total_TF_mRNA')
plt.title("colored by cell state")
plt.show()