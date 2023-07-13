#! /bin/bash

# This is a script to run convergence analysis on Laglos using Release version
laglos_dir="~/Workspace/Laglos"
create_convergence="~/Workspace/Laglos/scripts/create_convergence_table.py"
temp_output="~/Workspace/Laglos/saved/convergence/temp_output/"
convergence_dir="~/Workspace/Laglos/saved/convergence/"
scripts_dir="~/Workspace/Laglos/scripts"

cd $convergence_dir
rm *.out
rm *.txt

cd $temp_output
rm *.out
rm *.txt

cd $laglos_dir

# options="-m data/ref-tube.mesh -visc -mm -tf 0.225 -ot -cfl 0.01" # Shocktube testing, no rotation
options="-m data/ref-segment-tube.mesh -tf 0.225 -cfl 0.1 -ot -visc -mm -vis" # Shocktube testing, 1D
# options="-m ../data/rectangle_saltzman.mesh -visc -mm -tf 0.6 -ot -cfl 0.01" # Saltzman
# options="-m ../data/ref-square-c0.mesh -tf 2 -cfl 0.5 -ot -visc -mm" # Isentropic Vortex, moving center
# options="-m ../data/ref-square-c0.mesh -tf 0.2 -cfl 0.1 -ot -visc -mm"
# options="-m ../data/shocktube.mesh -tf .67 -cfl 0.5 -ot -visc -mm"
# Testing low order convergence isentropic vortex
# ./Laglos -rs 0 $options
# ./Laglos -rs 1 $options
./Laglos -rs 2 $options
./Laglos -rs 3 $options
./Laglos -rs 4 $options
./Laglos -rs 5 $options
./Laglos -rs 6 $options
./Laglos -rs 7 $options
./Laglos -rs 8 $options
# ./Laglos -rs 9 $options
# ./Laglos -rs 10 $options
# ./Laglos -rs 11 $options
cd $scripts_dir

# python stuff for iterating over files to build convergence table
python3 $create_convergence $temp_output
