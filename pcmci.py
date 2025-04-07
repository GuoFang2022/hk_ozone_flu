# %%
import pandas as pd
import os
import pyarrow.feather as feather
from matplotlib import pyplot as plt
from tigramite import data_processing as pp
from tigramite import plotting as tp
from tigramite.pcmci import PCMCI
from tigramite.independence_tests.parcorr import ParCorr
import numpy as np

# %%
data_all = pd.read_csv("df_pcmci_norm.csv")

# %%
data_all.rename(columns={'logitFluP': r'$Flu$',
                         'o3h8max': r'$O_{3}$',
                         'season': r'$Season$',
                         'temp': r'$T$',
                         'ah': r'$AH$'}, inplace=True)

# %%
data_all

# %%
# Preparing functions for specifying data format (mask or not) and variables under study ('var_names')
def get_data(data_all,  var_names=[]):

    data_out = data_all

    # Select the columns
    data_out = data_out[var_names]

    # Turn into numpy array
    data_out = data_out.values

    # Return
    return data_out

# %%
# Specify link assumptions
def define_link_assumptions(var_names, tau_min, tau_max):
    # Get index of the specified variable, if it exists
    if r'$T$' in var_names:
        temp_idx = np.argwhere(np.array(var_names) == r'$T$')[0, 0]
    else:
        temp_idx = None

    if r'$AH$' in var_names:
        ah_idx = np.argwhere(np.array(var_names) == r'$AH$')[0, 0]
    else:
        ah_idx = None
    
    # Build dictionary
    link_assumptions = {}
    
    for j, target_var in enumerate(var_names):
        # AH cannot be affected by flu and ozone at any lags
        if target_var == r'$AH$': 
            link_assumptions[j] = {(i, -tau):'-?>' for i, parent_var in enumerate(var_names)
                                   for tau in range(tau_min, tau_max+1) if 
                                   parent_var != r'$Flu$' and 
                                   parent_var != r'$O_{3}$' and 
                                   (i, -tau) != (j, 0)} 
        # T cannot be affected by flu and ozone at any lags, but can be affected by AH at non-zero lags
        elif target_var == r'$T$':
            link_assumptions[j] = {(i, -tau):'-?>' for i, parent_var in enumerate(var_names)
                                   for tau in range(tau_min, tau_max+1) if 
                                   parent_var != r'$Flu$' and 
                                   parent_var != r'$O_{3}$' and 
                                   (i, -tau) != (j, 0) and (i, -tau) != (ah_idx, 0)}
           
        # Ozone may be influenced by all variables except Flu
        elif target_var == r'$O_{3}$':
            link_assumptions[j] = {(i, -tau):'-?>' for i, parent_var in enumerate(var_names)
                                   for tau in range(tau_min, tau_max+1) if 
                                   parent_var != r'$Flu$' and (i, -tau) != (j, 0)}                    
                        
        # Default case for other variables (i.e., Flu, no restrictions, may be influenced by all variables at all lags)  
        else:
            link_assumptions[j] = {(i, -tau):'-?>' for i in range(len(var_names))
                                   for tau in range(tau_min, tau_max+1) if
                                    (i, -tau) != (j, 0)}
            
 
    return link_assumptions

# %%
def apply_pcmci(data_all,
                var_names,
                tau_min,
                tau_max,
                pc_alpha,
                verbosity):
    # Get the data 
    data = get_data(data_all=data_all, var_names=var_names)

    # Prepare the DataFrame object
    dataframe = pp.DataFrame(data, 
                             var_names=var_names,
                             missing_flag=999.)
    
    pcmci = PCMCI(dataframe=dataframe,
                  cond_ind_test=ParCorr(),
                  verbosity=verbosity)
    
    # Get the selected_links arguement
    link_assumptions =  define_link_assumptions(var_names, tau_min, tau_max)

    # Run PCMCI^+ with these parameters
    results = pcmci.run_pcmciplus(tau_min=tau_min,
                                  tau_max=tau_max,
                                  pc_alpha=pc_alpha, 
                                  link_assumptions= link_assumptions)
    
    
    # Plot summary causal graph
    tp.plot_graph(
        arrow_linewidth=7.0,
        figsize=(10, 5),
        vmin_edges=-0.4,
        vmax_edges=0.4,
        node_label_size=21,
        link_label_fontsize=18,
        val_matrix=results['val_matrix'],
        graph=results['graph'],
        var_names=var_names,
        link_colorbar_label='cross-MCI (edges)',
        node_colorbar_label='auto-MCI (nodes)',
        vmin_nodes=0,
        vmax_nodes=1.0, 
        cmap_nodes="Reds",
        node_aspect=None,
        node_size=0.3,
        tick_label_size= 12,
        label_fontsize=15,
        show_colorbar=1
    );
    
    # Automatically adjust the layout
    plt.tight_layout()
    plt.savefig(f"hk_flu_PCMCI_summary_graph_{pc_alpha}.jpg", format="jpg", dpi=300)
    plt.show()

    pcmci.print_significant_links(
        p_matrix = results['p_matrix'],
        val_matrix = results['val_matrix'],
        alpha_level = 0.05)

    # Print optimal alpha (if applicable)
    print(f"pc_alpha: {pc_alpha}")
    if "optimal_alpha" in results:
        print(f"Optimal alpha: {results['optimal_alpha']}")
  
    return results

# %%
# Set the parameters for the PCMCI analysis
pc_alpha= 0.05

results = apply_pcmci(data_all=data_all,
                      var_names=[r'$Flu$', r'$O_{3}$',r'$AH$', r'$T$'],
                      tau_min=0,
                      tau_max=2,
                      pc_alpha=pc_alpha,
                      verbosity=0
                     )

# Close the plot to free up memory
plt.close()

# %%
pc_alpha= 0.1

results = apply_pcmci(data_all=data_all,
                      var_names=[r'$Flu$', r'$O_{3}$',r'$AH$', r'$T$'],
                      tau_min=0,
                      tau_max=2,
                      pc_alpha=pc_alpha,
                      verbosity=0
                     )

# Close the plot to free up memory
plt.close()

# %%
pc_alpha= 0.01

results = apply_pcmci(data_all=data_all,
                      var_names=[r'$Flu$', r'$O_{3}$',r'$AH$', r'$T$'],
                      tau_min=0,
                      tau_max=2,
                      pc_alpha=pc_alpha,
                      verbosity=0
                     )

# Close the plot to free up memory
plt.close()


