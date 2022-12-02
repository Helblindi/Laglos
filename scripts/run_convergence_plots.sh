#! /bin/bash

# This is a script to run convergence analysis on ex9p-continuous
mfem_dir="/Users/sheridan7/Workspace/mfem/examples"
create_convergence_plot="/Users/sheridan7/Workspace/mfem/examples/ex9p-analysis/scripts/create_convergence_plot.py"
temp_output="/Users/sheridan7/Workspace/mfem/examples/ex9p-analysis/temp_output/"
analysis_dir="/Users/sheridan7/Workspace/mfem/examples/ex9p-analysis/"

cd $analysis_dir
rm *.out
rm *.txt

cd $temp_output
rm *.out
rm *.txt

cd $mfem_dir

# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 1 -rp 0 -p 4 -tf 2
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 2 -rp 0 -p 4 -tf 2
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 3 -rp 0 -p 4 -tf 2
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 4 -rp 0 -p 4 -tf 2
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 5 -rp 0 -p 4 -tf 2
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 6 -rp 0 -p 4 -tf 2
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 7 -rp 0 -p 4 -tf 2
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 8 -rp 0 -p 4 -tf 2


# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 1 -rp 0 -p 4 -ots
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 2 -rp 0 -p 4 -ots
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 3 -rp 0 -p 4 -ots
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 4 -rp 0 -p 4 -ots
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 5 -rp 0 -p 4 -ots
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 6 -rp 0 -p 4 -ots
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 7 -rp 0 -p 4 -ots
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 8 -rp 0 -p 4 -ots

# mpirun -np 4 ex9p-continuous -m ../data/periodic-square.mesh -dt .001 -rs 1 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/periodic-square.mesh -dt .001 -rs 2 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/periodic-square.mesh -dt .001 -rs 3 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/periodic-square.mesh -dt .001 -rs 4 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/periodic-square.mesh -dt .001 -rs 5 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/periodic-square.mesh -dt .001 -rs 6 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/periodic-square.mesh -dt .001 -rs 7 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/periodic-square.mesh -dt .001 -rs 8 -rp 0 -p 4 -tf 2.84

# mpirun -np 2 ex9p-continuous -m ../data/periodic-square.mesh -ct -rs 1 -rp 0 -p 4 -ots
# mpirun -np 2 ex9p-continuous -m ../data/periodic-square.mesh -ct -rs 2 -rp 0 -p 4 -ots
# mpirun -np 2 ex9p-continuous -m ../data/periodic-square.mesh -ct -rs 3 -rp 0 -p 4 -ots
# mpirun -np 2 ex9p-continuous -m ../data/periodic-square.mesh -ct -rs 4 -rp 0 -p 4 -ots
# mpirun -np 2 ex9p-continuous -m ../data/periodic-square.mesh -ct -rs 5 -rp 0 -p 4 -ots
# mpirun -np 2 ex9p-continuous -m ../data/periodic-square.mesh -ct -rs 6 -rp 0 -p 4 -ots
# mpirun -np 2 ex9p-continuous -m ../data/periodic-square.mesh -ct -rs 7 -rp 0 -p 4 -ots
# mpirun -np 2 ex9p-continuous -m ../data/periodic-square.mesh -ct -rs 8 -rp 0 -p 4 -ots

# mpirun -np 4 ex9p-continuous -m ../data/periodic-cube.mesh -ct -rs 1 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/periodic-cube.mesh -ct -rs 2 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/periodic-cube.mesh -ct -rs 3 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/periodic-cube.mesh -ct -rs 4 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/periodic-cube.mesh -ct -rs 5 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/periodic-cube.mesh -ct -rs 6 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/periodic-cube.mesh -ct -rs 7 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/periodic-cube.mesh -ct -rs 8 -rp 0 -p 4 -tf 2.84

cd $analysis_dir
mv *.out temp_output/.
cd scripts

# python stuff for iterating over files to build convergence table
python3 $create_convergence_plot $temp_output
