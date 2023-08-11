"""
This python script file plots output from multiple refinements.

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
num_cols = int(sys.argv[2])

def main():
   fig, axs = plt.subplots(1, num_cols)
   for filename in sorted(os.listdir(directory)):
      f = os.path.join(directory, filename)
      refinements = filename[3:5]
      print('refinements: ', refinements)
      # rec_arr = np.genfromtxt(f, delimiter = ",", names=True, dtype=None)
      df = pd.read_csv(f, dtype=float).sort_values(by=['x'])
      # print(df)
      # print(df.dtypes)
      # df = df.sort_values(by=['x'])
      # print('sorted df')
      # print(df)
      for i, col in enumerate(df.columns):
         if (i == 0): continue # Skip first column
         axs[i-1].set_xlabel(df.columns[0])
         axs[i-1].set_ylabel(col)
         axs[i-1].plot(df[df.columns[0]], df[col], label=refinements)
   
   plt.legend()
   plt.show()

# then we put main at the bottom to run everything
main()
        