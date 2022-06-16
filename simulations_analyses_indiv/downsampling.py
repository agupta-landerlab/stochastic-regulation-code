from sim_dependencies.py import *

## read in downsampling files, each corresponding to corrs for varying:
# capture efficiencies
# number of cells
# simulation run # for each ^ combo

capture_efficiency = ['10', '20', '30', '40', '50', '60', '70', '80', '90', '100'] #rows
num_cells = ['50k', '25k', '10k', '5k', '2.5k', '1k'] #columns

plt.figure(figsize=(20,15))

i=1
num_rows, num_cols = len(num_cells), len(capture_efficiency)
for NC_i in range(num_rows):
    for CE_i in range(num_cols):
        NC = num_cells[NC_i]
        CE = capture_efficiency[-num_cols:][CE_i]
        
        pair_df = pd.DataFrame()
        simruns=[]
        print('one more CE+CN pair')
        for simrun in range(0,25):
            try:
                fn = 'simulated_data_for_figures/fig4C_downsampling/for_metaplots/CE='+CE+'_NC='+NC+"_corrs_run"+str(simrun)+'.csv.csv'
                these_corrs = pd.DataFrame(pd.read_csv(fn))
                simruns.extend([simrun]*len(these_corrs))
                pair_df = pd.concat([pair_df, these_corrs], axis=0)
                pair_df['sim_run'] = simruns
            except:
                pass
        
        #DEFINING # SUBPLOTS: (1) NUM ROWS (2) NUM COLS (3) INDEX FOR THIS SUBPLOT
        #SOURCE: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplot.html
        ax = plt.subplot(num_rows, num_cols, i)#, frameon=False)
    
        try:
            sns.set_style("whitegrid")
            sns.stripplot(data=pair_df,x='timepoint',y='TF_total_t0:Target_nascent_tN',
                          color='black',size=2,jitter=False)
            plt.plot(np.arange(0,48,1), list(pair_df.groupby("timepoint").mean()['TF_total_t0:Target_nascent_tN']),
                     color='white')
            plt.xticks(" ")
            plt.yticks(" ")
            plt.ylabel(" ")
            plt.xlabel(" ")
            sns.despine()
    
            #change spacing between subplots
            plt.subplots_adjust(wspace=0.05, hspace=0.1)
        except:
            pass
        
        i+=1

plt.show()

## for main figure (fewer values plotted)

capture_efficiency = ['10', '50','90'] #rows
num_cells = ['25k', '5k', '1k'] #columns

fig_main_DS = plt.figure(figsize=(20,8))

i=1
num_rows, num_cols = len(num_cells), len(capture_efficiency)
for NC_i in range(num_rows):
    for CE_i in range(num_cols):
        NC = num_cells[NC_i]
        CE = capture_efficiency[-num_cols:][CE_i]
        
        pair_df = pd.DataFrame()
        simruns=[]
        print('one more CE+CN pair')
        for simrun in range(0,25):
            try:
                fn = 'simulated_data_for_figures/fig4C_downsampling/for_metaplots/CE='+CE+'_NC='+NC+"_corrs_run"+str(simrun)+'.csv.csv'
                these_corrs = pd.DataFrame(pd.read_csv(fn))
                simruns.extend([simrun]*len(these_corrs))
                pair_df = pd.concat([pair_df, these_corrs], axis=0)
                pair_df['sim_run'] = simruns
            except:
                pass
        
        #DEFINING # SUBPLOTS: (1) NUM ROWS (2) NUM COLS (3) INDEX FOR THIS SUBPLOT
        #SOURCE: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplot.html
        ax = plt.subplot(num_rows, num_cols, i)#, frameon=False)
    
        try:
            sns.set_style("whitegrid")
            sns.stripplot(data=pair_df,x='timepoint',y='TF_total_t0:Target_nascent_tN',
                          color='#6f63ac',size=3,jitter=False)
            plt.plot(np.arange(0,48,1), list(pair_df.groupby("timepoint").mean()['TF_total_t0:Target_nascent_tN']),
                     color='#6f63ac', linewidth=4)
            plt.yticks(" ")
            plt.ylabel(" ")
            plt.xlabel(" ")
            sns.despine()
            #change spacing between subplots
            plt.subplots_adjust(wspace=0.05, hspace=0.1)
        except:
            pass
        
        i+=1