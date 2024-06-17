#! /bin/bash

# This is a script to run convergence analysis on Laglos using Release version
# These values are specific to Whistler
laglos_dir="/home/sheridanm/software/Laglos"
scripts_dir="${laglos_dir}/scripts"
bin_dir="${laglos_dir}/build"
create_convergence="${scripts_dir}/create_convergence_table.py"
results_dir="/home/sheridanm/scratch/Laglos-results"
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

# EDIT THESE PARAMS
mesh_file="data/ref-square-c0.mesh"
problem=4
final_time=.6
cfl=0.5
mv_option=2
fv_option=2
mv_it_op=2
mv_iter_n=2
output_location="Notes/cellvol-20240614/it2/noh/cflp5-novisc"
output_file="${bin_dir}/out-noh-novisc-r"
########

options="-m ${mesh_file} -p ${problem} "
options+="-tf ${final_time} -cfl ${cfl} "
options+="-mv ${mv_option} -fv ${fv_option} "
options+="-mv-it-op ${mv_it_op} -mv-iter-n ${mv_iter_n} " 
options+="-of ${output_location} "

echo $options

#srun-serial -n 1 ./Laglos -rs 1 $options > "${output_file}1" &
#srun-serial -n 1 ./Laglos -rs 2 $options > "${output_file}2" &
srun-serial -n 1 ./Laglos -rs 3 $options > "${output_file}3" &
srun-serial -n 1 ./Laglos -rs 4 $options > "${output_file}4" &
srun-serial -n 1 ./Laglos -rs 5 $options > "${output_file}5" -print &
srun-serial -n 1 ./Laglos -rs 6 $options > "${output_file}6" &
srun-serial -n 1 ./Laglos -rs 7 $options > "${output_file}7" &
srun-serial -n 1 ./Laglos -rs 8 $options > "${output_file}8" &

# Lastly, copy this script to the output_file location
# for documentation on run parameters
cp "${scripts_dir}/whistler_convergence_script.sh" "${results_dir}/${output_location}/." 

# cd $scripts_dir

# python stuff for iterating over files to build convergence table
# python3 $create_convergence $temp_output

