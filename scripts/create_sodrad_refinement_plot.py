"""
This python script file plots output from multiple refinements.

The second command line argument corresponds to which quantity 
the user would like to plot.  These correspond to the following:
   1: Density
   2: Velocity in x direction
   3: Specific Total Energy
   4: Pressure
   5: Sound speed

Note: The label in the plot must be modified according to the mesh 
      size of the problem.

Example:
> python3 create_refinement_plot.py "/Users/madisonsheridan/Workspace/Laglos/build/results/state_vectors/"
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
assert len(sys.argv) == 2, "This file needs 1 input arguments: directory, but " + str(len(sys.argv)) + " were given."
directory = str(sys.argv[1])
screen_aspect = 0.75 # y_length / x_length

def main():
   fig, ax = plt.subplots(figsize=(13.5,10)) # Gives ~0.75 aspect ratio
   ###############
   ### Density ###
   ###############
   for filename in sorted(os.listdir(directory)):
      # if flag % 1 == 0:
      f = os.path.join(directory, filename)
      refinements = filename[3:5]
      
      df = pd.read_csv(f).sort_values(by=['x'])
      df = df.drop(columns=['cell_type'], axis=1)

      if (refinements != "ex"):
         _label = "# dofs = " + str(df.shape[0])
         ax.scatter(df[df.columns[0]], df[df.columns[1]], label=_label, s=2)
      

   # Exact solution is stored in a file in data directory
   # Exact solution columns are x, rho, u, p
   sod_rad_file = "../data/sodrad-sol-EP2D-Osher-00589.txt"
   df_exact = pd.read_csv(sod_rad_file, delim_whitespace=True).sort_values(by=['x'])
   ax.plot(df_exact[df_exact.columns[0]], df_exact[df_exact.columns[1]], 'k-', label="Reference solution")
   
   # ax.gca().set_aspect(screen_aspect / ax.gca().get_data_ratio())
   ax.legend(markerscale=4)
   ax.set_xlim(0., 1.)
   ax.set_xlabel("$x$")
   ax.set_ylabel("Density")
   _dpi = fig.get_dpi()
   fig.savefig("sodrad-density.png", dpi=_dpi)

   #########################
   ##### Velocity norm #####
   #########################
   ax.clear()
   for filename in sorted(os.listdir(directory)):
      # if flag % 1 == 0:
      f = os.path.join(directory, filename)
      refinements = filename[3:5]
      
      df = pd.read_csv(f).sort_values(by=['x'])
      df = df.drop(columns=['cell_type'], axis=1)

      if (refinements != "ex"):
         _label = "# dofs = " + str(df.shape[0])
         ax.scatter(df[df.columns[0]], df[df.columns[2]], label=_label, s=2)
      

   # Exact solution is stored in a file in data directory
   # Exact solution columns are x, rho, u, p
   sod_rad_file = "../data/sodrad-sol-EP2D-Osher-00589.txt"
   df_exact = pd.read_csv(sod_rad_file, delim_whitespace=True).sort_values(by=['x'])
   ax.plot(df_exact[df_exact.columns[0]], df_exact[df_exact.columns[2]], 'k-', label="Reference solution")
   
   # ax.gca().set_aspect(screen_aspect / ax.gca().get_data_ratio())
   ax.legend(markerscale=4)
   ax.set_xlim(0., 1.)
   ax.set_xlabel("$x$")
   ax.set_ylabel("$||v||$")
   _dpi = fig.get_dpi()
   fig.savefig("sodrad-velocity.png", dpi=_dpi)

   ####################
   ##### Pressure #####
   ####################
   ax.clear()
   for filename in sorted(os.listdir(directory)):
      # if flag % 1 == 0:
      f = os.path.join(directory, filename)
      refinements = filename[3:5]
      
      df = pd.read_csv(f).sort_values(by=['x'])
      df = df.drop(columns=['cell_type'], axis=1)

      if (refinements != "ex"):
         _label = "# dofs = " + str(df.shape[0])
         ax.scatter(df[df.columns[0]], df[df.columns[4]], label=_label, s=2)
      

   # Exact solution is stored in a file in data directory
   # Exact solution columns are x, rho, u, p
   sod_rad_file = "../data/sodrad-sol-EP2D-Osher-00589.txt"
   df_exact = pd.read_csv(sod_rad_file, delim_whitespace=True).sort_values(by=['x'])
   ax.plot(df_exact[df_exact.columns[0]], df_exact[df_exact.columns[3]], 'k-', label="Reference solution")
   
   # ax.gca().set_aspect(screen_aspect / ax.gca().get_data_ratio())
   ax.legend(markerscale=4)
   ax.set_xlim(0., 1.)
   ax.set_xlabel("$x$")
   ax.set_ylabel("Pressure")
   _dpi = fig.get_dpi()
   fig.savefig("sodrad-pressure.png", dpi=_dpi)

# then we put main at the bottom to run everything
main()
        