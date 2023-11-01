"""
This python script is a workaround for the lack of any sort of for loop
in the glvis script syntax. This file will help generate a GLVis 
animation script when there are far too many output files to tediously
add a screenshot line for each.

Example:
> python create_glvis_animation_script.py "/Users/madisonsheridan/Workspace/Laglos/build/vortex/gfs_r03/"
"""

import numpy as np
import re
import os
import sys


# check command line inputs
assert len(sys.argv) == 2, "This file needs 1 input argument: directory, but " + str(len(sys.argv)) + " were given."
directory_path = str(sys.argv[1])
filenames_without_suffix = []

# now we define the main function to be called at the end
def main():
   # comment out "run_simuations" if you only want to compute the errors
   vals = gather_files()
   write_glvis_script("rho")
   write_glvis_script("v")
   write_glvis_script("ste")


def gather_files():
   # Check if the directory exists
   if os.path.exists(directory_path) and os.path.isdir(directory_path):
      # List all files in the directory
      files = sorted(os.listdir(directory_path))
      # Loop through the files
      for filename in files:
         # Check if the file has a .txt suffix
         if filename.endswith('.mesh'):
               # Remove the .mesh suffix and append to the list
               filenames_without_suffix.append(filename[:-5])


def write_glvis_script(_flag):
   print("writing glvis script for ", _flag)
   glvis_filename = directory_path + _flag + ".glvis"
   f = open(glvis_filename, "w")

   f.write("# Visualization window geometry\n")
   f.write("window 0 0 900 900\n")
   f.write("\n# Initial solution\n")

   line_str = "solution " + filenames_without_suffix[0] + ".mesh " + filenames_without_suffix[0] + "_" + _flag + ".gf\n\n"
   f.write(line_str)
   f.write("# Setup the GLVis scene. Executed after pressing the space bar.\n")
   f.write("{\n")
   f.write("\tperspective off\n")
   f.write("\tview 0 0\n")
   f.write("\tviewcenter 0 0\n")
   f.write("\tvaluerange 0.0 1.0\n")
   f.write("\tautoscale off\n")
   f.write("\tzoom 1.8\n")
   f.write("\tkeys rRmclppppppppppppppppppppppppppp\n")
   f.write("}\n")
   f.write("# Take multiple screenshots. Executed after pressing the space bar.\n")
   f.write("{\n")
   for filename in filenames_without_suffix:
       line_str = "\tsolution " + filename + ".mesh " + filename + "_" + _flag + ".gf screenshot " + filename + "_" + _flag + ".png\n"
       f.write(line_str)
   f.write("}\n")

   f.close()

# then we put main at the bottom to run everything
main()
