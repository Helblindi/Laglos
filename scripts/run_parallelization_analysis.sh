#! /bin/bash

# This is a script to run convergence analysis on ex9p-continuous
mfem_dir="/Users/madisonsheridan/Workspace/mfem/examples"
create_parallelization="/Users/madisonsheridan/Workspace/mfem/examples/ex9p-analysis/scripts/create_parallelization_table.py"
temp_output="/Users/madisonsheridan/Workspace/mfem/examples/ex9p-analysis/temp_output/"
analysis_dir="/Users/madisonsheridan/Workspace/mfem/examples/ex9p-analysis/"
cd $analysis_dir
rm *.out
rm *.txt

cd $temp_output
rm *.out
rm *.txt

cd $mfem_dir
# mpirun -np 1 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 7 -rp 0 -p 4 -tf 2
# mpirun -np 2 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 7 -rp 0 -p 4 -tf 2
# mpirun -np 4 ex9p-continuous -m ../data/periodic-segment.mesh -ct -rs 7 -rp 0 -p 4 -tf 2

# mpirun -np 1 ex9p-continuous -m ../data/ref-square.mesh -ct -rs 4 -rp 0 -p 4 -tf 2.84
# mpirun -np 2 ex9p-continuous -m ../data/ref-square.mesh -ct -rs 4 -rp 0 -p 4 -tf 2.84
# mpirun -np 4 ex9p-continuous -m ../data/ref-square.mesh -ct -rs 4 -rp 0 -p 4 -tf 2.84

mpirun -np 1 ex9p-continuous -m ../data/periodic-square.mesh -ot -rs 7 -rp 0 -p 4 -tf 2.83
mpirun -np 2 ex9p-continuous -m ../data/periodic-square.mesh -ot -rs 7 -rp 0 -p 4 -tf 2.83
mpirun -np 4 ex9p-continuous -m ../data/periodic-square.mesh -ot -rs 7 -rp 0 -p 4 -tf 2.83

cd $analysis_dir
mv *.out temp_output/.
cd scripts

# python stuff for iterating over files to build convergence table
python3 $create_parallelization $temp_output
