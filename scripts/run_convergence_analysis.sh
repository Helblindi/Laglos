#! /bin/bash

# This is a script to run convergence analysis on Laglos using Release version
laglos_dir="/Users/madisonsheridan/Workspace/Laglos"
scripts_dir="${laglos_dir}/scripts"
bin_dir="${laglos_dir}/build"
create_convergence="${scripts_dir}/create_convergence_table.py"
create_refinement="${scripts_dir}/create_refinement_plot.py"
temp_output="${bin_dir}/results/convergence/temp_output/"
state_vectors="${bin_dir}/results/state_vectors/"
convergence_dir="${bin_dir}/results/convergence/"


# cd $convergence_dir
# rm *.out
# rm *.txt

# cd $temp_output
# rm -rf *
# cd $state_vectors
# rm -rf *

cd $bin_dir

# Dim = 1 cases
# options="-m data/ref-segment.mesh -visc -mm -tf 0.6 -ot -cfl 1 -so -vs 100" # 1D Smooth Wave 2nd Order IDP Paper
# options="-m data/ref-segment.mesh -p 1 -tf 0.225 -cfl 0.25 -ot -visc -mm -so -vs 100" # Sod
# options="-m data/ref-segment.mesh -p 2 -tf 0.15 -cfl 0.5 -ot -visc -mm -so -vs 100"   ## Lax
# options="-m data/ref-segment.mesh -p 3 -tf 0.667 -cfl 0.2 -ot -visc -mm -so -vs 100"    ## Leblanc


# Dim = 2 cases
options="-m data/square5c0_vortex.mesh -p 5 -tf 4 -cfl 0.01 -ot -visc -mm -so -vs 100 -of isen-vortex" ## Isentropic Vortex
# options="-m data/distorted-square.mesh -p 1 -tf 0.225 -cfl 0.01 -ot -visc -mm -so -vs 100 -of sod-distorted" #rs0-4 ## Sod on distorted grid
# options="-m data/ref-square.mesh -p 1 -tf 0.225 -cfl 0.5 -ot -visc -mm -so -vs 100 -of sod-cartesian" #rs2-6 ## Sod on cartesian grid

./Laglos -rs 0 $options
./Laglos -rs 1 $options
./Laglos -rs 2 $options
./Laglos -rs 3 $options
./Laglos -rs 4 $options
# ./Laglos -rs 5 $options
# ./Laglos -rs 6 $options
# ./Laglos -rs 7 $options
# ./Laglos -rs 8 $options
# ./Laglos -rs 9 $options
# ./Laglos -rs 10 $options
# ./Laglos -rs 11 $options
# ./Laglos -rs 12 $options

# cd $scripts_dir

# python stuff for iterating over files to build convergence table
# python3 $create_convergence $temp_output
# python3 $create_refinement $state_vectors 2
