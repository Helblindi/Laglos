#! /bin/bash

# This is a script to run convergence analysis on Laglos using Release version
laglos_dir="/Users/madisonsheridan/Workspace/Laglos"
scripts_dir="${laglos_dir}/scripts"
bin_dir="${laglos_dir}/build"
create_convergence="${scripts_dir}/create_convergence_table.py"
temp_output="${bin_dir}/results/convergence/temp_output/"
convergence_dir="${bin_dir}/results/convergence/"


cd $convergence_dir
rm *.out
rm *.txt

cd $temp_output
rm -rf *

cd $bin_dir

# Dim = 1 cases
# options="-m data/ref-segment.mesh -visc -mm -tf 0.6 -ot -cfl 1 -so -vs 100" # 1D Smooth Wave 2nd Order IDP Paper
# options="-m data/ref-segment.mesh -visc -mm -tf 0.225 -ot -cfl 0.25 -so -vs 100" # Shocktube testing, no rotation

# Dim = 2 cases
options="-m data/square5c0_vortex.mesh -tf 2 -cfl 0.5 -ot -visc -mm -vs 100 -so" # Isentropic Vortex, stationary center

# ./Laglos -rs 0 $options
./Laglos -rs 1 $options
./Laglos -rs 2 $options
./Laglos -rs 3 $options
./Laglos -rs 4 $options
./Laglos -rs 5 $options
# ./Laglos -rs 6 $options
# ./Laglos -rs 7 $options
# ./Laglos -rs 8 $options
# ./Laglos -rs 9 $options
# ./Laglos -rs 10 $options
# ./Laglos -rs 11 $options
# ./Laglos -rs 12 $options
cd $scripts_dir

# python stuff for iterating over files to build convergence table
python3 $create_convergence $temp_output
