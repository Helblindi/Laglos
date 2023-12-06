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
> python3 create_refinement_plot.py "/Users/madisonsheridan/Workspace/Laglos/build/results/state_vectors/ 2"
"""

import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
from tabulate import tabulate


# check command line inputs
assert len(sys.argv) == 3, "This file needs 2 input arguments: directory, num cols, but " + str(len(sys.argv)) + " were given."
directory = str(sys.argv[1])
col = int(sys.argv[2])

def main():
   # fig, axs = plt.subplots(1, 1)
   flag = 1
   for filename in sorted(os.listdir(directory)):
      if flag % 2 == 0:
         f = os.path.join(directory, filename)
         refinements = filename[3:5]
         
         print('refinements: ', refinements)
         df = pd.read_csv(f, dtype=float).sort_values(by=['x'])

         # if (i == 0): continue # Skip first column
         # plt.set_xlabel(df.columns[0])
         # plt.set_ylabel("$\rho$")
         if (refinements == "ex"):
            plt.plot(df[df.columns[0]], df[df.columns[1]], 'k-', label="Exact solution")
         else:
            num_ref = int(refinements)
            _label = "# dof = " + str(1 * pow(2,num_ref))
            plt.plot(df[df.columns[0]], df[df.columns[col]], '-', label=_label, linewidth=0.5)

         # plt.xlabel(df.columns[0])
         # plt.ylabel(df.columns[1])
         # plt.title("Sod Shocktube")

         # for i, col in enumerate(df.columns):
         #    if (i == 0): continue # Skip first column
         #    axs[i-1].set_xlabel(df.columns[0])
         #    axs[i-1].set_ylabel(col)
         #    if (refinements == "ex"):
         #       axs[i-1].plot(df[df.columns[0]], df[col], 'k--', label=refinements)
         #    else:
         #       axs[i-1].plot(df[df.columns[0]], df[col], '-', label=refinements)

         # plt.xlabel(df.columns[0])
         # plt.ylabel(df.columns[1])
         # plt.title("Sod Shocktube")
      
      flag += 1
   
   plt.legend()
   plt.show()

# then we put main at the bottom to run everything
main()
        