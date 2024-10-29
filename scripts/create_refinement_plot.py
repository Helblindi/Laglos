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
   flag = 0
   for filename in sorted(os.listdir(directory)):
      # if flag % 1 == 0:
      f = os.path.join(directory, filename)
      refinements = filename[3:5]
      
      print('refinements: ', refinements)
      df = pd.read_csv(f).sort_values(by=['x'])
      df = df.drop(columns=['cell_type'], axis=1)

      # if (i == 0): continue # Skip first column
      # plt.set_xlabel(df.columns[0])
      # plt.set_ylabel("$\rho$")
      if (refinements == "ex"):
         flag -= 1
         plt.plot(df[df.columns[0]], df[df.columns[col]], 'k-', label="Analytical solution")
      else:
         # num_ref = int(refinements)
         # _label = str(1 * pow(2,num_ref))
         # _label = refinements
         _label = "# dof = " + str(df.shape[0])
         plt.scatter(df[df.columns[0]], df[df.columns[col]], label=_label, s=0.3)
      
      # flag += 1
      
      
   
   screen_aspect = 0.75 # y_length / x_length
   plt.gca().set_aspect(screen_aspect / plt.gca().get_data_ratio())
   plt.legend(markerscale=4)
   plt.xlim(0., 1.1)
   plt.show()

# then we put main at the bottom to run everything
main()
        