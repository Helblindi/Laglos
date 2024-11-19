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
assert len(sys.argv) == 3, "This file needs 2 input arguments: directory, num cols, but " + str(len(sys.argv)) + " were given."
directory = str(sys.argv[1])
col = int(sys.argv[2])

def main():
   fig, ax = plt.subplots(figsize=(10,7)) # Gives ~0.75 aspect ratio
   flag = 0
   for filename in sorted(os.listdir(directory)):
      if (filename[-4:] == ".out"):
         f = os.path.join(directory, filename)
         indic = filename[0:-4]
         print("indic: ", indic)
         df = pd.read_csv(f).sort_values(by=['x'])
         df = df.drop(columns=['cell_type'], axis=1)

         if (indic == "0ppd"):
            ax.scatter(df[df.columns[0]], df[df.columns[col]], label="Mass conservative", s=2)   
         elif (indic == "1mv2"):
            ax.scatter(df[df.columns[0]], df[df.columns[col]], label="Non conservative", s=2)
         elif (indic == "2exact"):
            ax.plot(df[df.columns[0]], df[df.columns[col]], 'k-', label="Exact solution") 
   
   ax.legend()
   ax.set_ylim(0., 18)
   ax.set_yticks(range(0, 19, 2))  # Start at 0, stop at 16, step by 2
   ax.set_xlim(0., 1)
   ax.set_xlabel("$x$")
   # ax.set_ylabel("$x$-Velocity")
   ax.set_ylabel("Density")
   _dpi = fig.get_dpi()
   fig.savefig("Fig1py.png", dpi=_dpi)
   plt.show()

# then we put main at the bottom to run everything
main()
        