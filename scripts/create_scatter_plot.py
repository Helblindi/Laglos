"""
This python script file plots output from a single refinement.  
It is mainly used to compare a single spatial refinement approximation
to the exact solution.  

Example:
> python3 create_scatter_plot.py /Users/madisonsheridan/Workspace/Laglos/build/results/state_vectors/sv_05.out
                                 /Users/madisonsheridan/Workspace/Laglos/build/results/state_vectors/sv_exact.out
                                 num_dofs
"""

import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

# Globally set parameters for figures
params = {'legend.markerscale': 10,
          'legend.fontsize': 24,
          'axes.labelsize': 24,
          'axes.titlesize': 24,
          'xtick.labelsize': 24,
          'ytick.labelsize': 24}
pylab.rcParams.update(params)

# check command line inputs
# assert len(sys.argv) == 4, "This file needs 3 input arguments: directory, num cols, but " + str(len(sys.argv)) + " were given."
approx_loc = str(sys.argv[1])
exact_loc = str(sys.argv[2])
dofs = int(sys.argv[3])

def main():
   fig, ax = plt.subplots(figsize=(13.5,10)) # Gives ~0.75 aspect ratio

   # Read in interior cells from dataframe and sort by x
   df_approx = pd.read_csv(approx_loc).query('cell_type != "bdr"').sort_values(by=['x'])
   df_exact = pd.read_csv(exact_loc).query('cell_type != "bdr"').sort_values(by=['x'])

   # Drop column corresponding to cell_type
   df_approx = df_approx.drop(columns=['cell_type'], axis=1)
   df_exact = df_exact.drop(columns=['cell_type'], axis=1)
   
   # Plot
   ax.plot(df_exact[df_exact.columns[0]], df_exact[df_exact.columns[1]], 'k-', label="Exact solution")
   _label = "# dof = " + str(dofs)
   ax.scatter(df_approx[df_approx.columns[0]], df_approx[df_approx.columns[1]], label=_label, s=.1)
   
   # ax.legend()
   # plt.show()
   dpi = fig.get_dpi()
   fig.savefig("Fig1py.png", dpi=dpi)

# then we put main at the bottom to run everything
main()
        