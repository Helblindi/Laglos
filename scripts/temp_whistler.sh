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


cd $bin_dir

# EDIT THESE PARAMS
mesh_file="data/square5c0_vortex.mesh"
problem=5
final_time=50
mv_option=2
fv_option=1
#mv_it_op=3
#mv_iter_n=3
output_location="Notes/mv-iter-analysis-20240503/time-of-collapse/"
output_file="${bin_dir}/out-ivtoc-"
########

options="-m ${mesh_file} -p ${problem} "
options+="-tf ${final_time} -rs 6 -print "
options+="-mv ${mv_option} -fv ${fv_option} "
#options+="-of ${output_location} "

echo $options

#it 0
#srun-serial ./Laglos $options -cfl 0.5 -of "${output_location}it0-cflp5" > "${output_file}it0-cflp5" &
#srun-serial ./Laglos $options -cfl 0.25 -of "${output_location}it0-cflp25" > "${output_file}it0-cflp25" &

# it 1
#srun-serial ./Laglos -mv-iter-n 1 -mv-it-op 1 -cfl 0.5 $options -of "${output_location}it1-op1-cflp5" > "${output_file}it1-op1-cflp5" &
#srun-serial ./Laglos -mv-iter-n 1 -mv-it-op 1 -cfl 0.25 $options -of "${output_location}it1-op1-cflp25" > "${output_file}it1-op1-cflp25" &

#srun-serial ./Laglos -mv-iter-n 1 -mv-it-op 2 -cfl 0.5 $options -of "${output_location}it1-op2-cflp5" > "${output_file}it1-op2-cflp5" &
#srun-serial ./Laglos -mv-iter-n 1 -mv-it-op 2 -cfl 0.25 $options -of "${output_location}it1-op2-cflp25" > "${output_file}it1-op2-cflp25" &

#srun-serial ./Laglos -mv-iter-n 1 -mv-it-op 3 -cfl 0.5 $options -of "${output_location}it1-op3-cflp5" > "${output_file}it1-op3-cflp5" &
#srun-serial ./Laglos -mv-iter-n 1 -mv-it-op 3 -cfl 0.25 $options -of "${output_location}it1-op3-cflp25" > "${output_file}it1-op3-cflp25" &

# it3
#srun-serial ./Laglos -mv-iter-n 3 -mv-it-op 1 -cfl 0.5 $options -of "${output_location}it3-op1-cflp5" > "${output_file}it3-op1-cflp5" &
#srun-serial ./Laglos -mv-iter-n 3 -mv-it-op 1 -cfl 0.25 $options -of "${output_location}it3-op1-cflp25" > "${output_file}it3-op1-cflp25" &

#srun-serial ./Laglos -mv-iter-n 3 -mv-it-op 2 -cfl 0.5 $options -of "${output_location}it3-op2-cflp5" > "${output_file}it3-op2-cflp5" &
#srun-serial ./Laglos -mv-iter-n 3 -mv-it-op 2 -cfl 0.25 $options -of "${output_location}it3-op2-cflp25" > "${output_file}it3-op2-cflp25" &

#srun-serial ./Laglos -mv-iter-n 3 -mv-it-op 3 -cfl 0.5 $options -of "${output_location}it3-op3-cflp5" > "${output_file}it3-op3-cflp5" &
#srun-serial ./Laglos -mv-iter-n 3 -mv-it-op 3 -cfl 0.25 $options -of "${output_location}it3-op3-cflp25" > "${output_file}it3-op3-cflp25" &

# it5
srun-serial ./Laglos -mv-iter-n 5 -mv-it-op 1 -cfl 0.5 $options -of "${output_location}it5-op1-cflp5" > "${output_file}it5-op1-cflp5" &
srun-serial ./Laglos -mv-iter-n 5 -mv-it-op 1 -cfl 0.25 $options -of "${output_location}it5-op1-cflp25" > "${output_file}it5-op1-cflp25" &

srun-serial ./Laglos -mv-iter-n 5 -mv-it-op 2 -cfl 0.5 $options -of "${output_location}it5-op2-cflp5" > "${output_file}it5-op2-cflp5" &
srun-serial ./Laglos -mv-iter-n 5 -mv-it-op 2 -cfl 0.25 $options -of "${output_location}it5-op2-cflp25" > "${output_file}it5-op2-cflp25" &

srun-serial ./Laglos -mv-iter-n 5 -mv-it-op 3 -cfl 0.5 $options -of "${output_location}it5-op3-cflp5" > "${output_file}it5-op3-cflp5" &
srun-serial ./Laglos -mv-iter-n 5 -mv-it-op 3 -cfl 0.25 $options -of "${output_location}it5-op3-cflp25" > "${output_file}it5-op3-cflp25" &

#it 10
srun-serial ./Laglos -mv-iter-n 10 -mv-it-op 1 -cfl 0.5 $options -of "${output_location}it10-op1-cflp5" > "${output_file}it10-op1-cflp5" &
srun-serial ./Laglos -mv-iter-n 10 -mv-it-op 1 -cfl 0.25 $options -of "${output_location}it10-op1-cflp25" > "${output_file}it10-op1-cflp25" &

srun-serial ./Laglos -mv-iter-n 10 -mv-it-op 2 -cfl 0.5 $options -of "${output_location}it10-op2-cflp5" > "${output_file}it10-op2-cflp5" &
srun-serial ./Laglos -mv-iter-n 10 -mv-it-op 2 -cfl 0.25 $options -of "${output_location}it10-op2-cflp25" > "${output_file}it10-op2-cflp25" &

srun-serial ./Laglos -mv-iter-n 10 -mv-it-op 3 -cfl 0.5 $options -of "${output_location}it10-op3-cflp5" > "${output_file}it10-op3-cflp5" &
srun-serial ./Laglos -mv-iter-n 10 -mv-it-op 3 -cfl 0.25 $options -of "${output_location}it10-op3-cflp25" > "${output_file}it10-op3-cflp25" &















