#plotting gene-specific kinetic parameter values for table 1 (literature-derived)
#interquartile ranges for each intrinsic parameter, taken from experiments reported in literature
#(refs cited in the main text)

#time between bursts (1/kon) (hours)
interburst_times = ([2.13, 3.7, 7.14], [2.4, 4, 7.7])

#burst duration (1/koff) (minutes)
burst_durations = ([1.25, 7.14, 39.6], [0.78, 7.8, 29.4])

#mRNA half-life (hours)
mRNA_thalfs = ([1.4, 2.5, 3.6], [2.67, 3.7, 5.16])

#protein half-life (hours)
prot_thalfs = ([15,28,62], [24,48,88])

#burst size
burst_sizes = ([3.2,5.2,9.6], [3.78,6.6,12.2])

#translation rate
translation_rates = ([22,59,179], [37,117,309])

all_intrinsic_params=[interburst_times,burst_durations,
                     mRNA_thalfs,prot_thalfs,burst_sizes,
                     translation_rates]

for TF_vals, nonTF_vals in all_intrinsic_params:
    TF_Q1, TF_med, TF_Q3 = TF_vals
    nonTF_Q1, nonTF_med, nonTF_Q3 = nonTF_vals
    lw = 4
    ms = 5
    fs = 20

    fig1 = plt.figure(figsize=(5,2))
    ax1 = plt.axes(frameon=False)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    plt.plot([TF_Q1,TF_Q1],[0,1],color='DarkBlue',linewidth=lw)
    plt.plot([TF_med,TF_med],[0,2],color='DarkBlue',linewidth=lw)
    plt.plot([TF_Q3,TF_Q3],[0,1],color='DarkBlue',linewidth=lw)
    plt.plot([nonTF_Q1,nonTF_Q1],[0,-1],color='Red',linewidth=lw)
    plt.plot([nonTF_med,nonTF_med],[0,-2],color='Red',linewidth=lw)
    plt.plot([nonTF_Q3,nonTF_Q3],[0,-1],color='Red',linewidth=lw)
    
    for param_val in [TF_Q1, TF_Q3]:
        plt.text(param_val,1.75,str(param_val),horizontalalignment='center',fontsize=fs,fontweight='bold')
        
    for param_val in [nonTF_Q1, nonTF_Q3]:
        plt.text(param_val,-1.75,str(param_val),horizontalalignment='center',fontsize=fs,fontweight='bold')
        
    plt.text(TF_med,2.75,str(TF_med),horizontalalignment='center',fontsize=fs,fontweight='bold')
    plt.text(nonTF_med,-2.75,str(nonTF_med),horizontalalignment='center',fontsize=fs,fontweight='bold')
    
    plt.plot([TF_Q1-.5,nonTF_Q3+1],[0,0],color='black',linewidth=2)
    sns.despine()
    plt.show()