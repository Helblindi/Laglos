#! /bin/bash

# This is a script to run all sample problems using a similar number of mesh cells on Laglos using Release version
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
#mv_option=2
#fv_option=1
#mv_it_op=3
#mv_iter_n=0
output_location="testing/mv-2-nocm-comparison-20240827/"
########

output_files_dir="${results_dir}/${output_location}/output_files/"
mkdir "${results_dir}/${output_location}"
mkdir "${output_files_dir}"

#options+=" -mv ${mv_option} -fv ${fv_option} "
#options+=" -print -of ${output_location} "
options+=" -mv 2 -fv 1 -no-cm " # Old method that results in parachuting
#options+=" -mv 8 -fv 1 -mv-iter-n 3 " #hiop
echo $options

# Sod
srun-serial ./laglos $options -of "${output_location}/sod" -m data/ref-square.mesh -p 1 -tf 0.225 -cfl 0.5 -rs 5 > "${output_files_dir}/out-sod" &

# Sod distorted
srun-serial ./laglos $options -of "${output_location}/sod-dis" -m data/distorted-square.mesh -p 1 -tf 0.225 -cfl 0.5 -rs 3 > "${output_files_dir}/out-sod-dis" &

# Sod radial
srun-serial ./laglos $options -of "${output_location}/sod-rad" -m data/ref-square.mesh -p 13 -tf 0.2 -cfl 0.25 -rs 5 > "${output_files_dir}/out-sod-rad" &

# Sod tube
# Sod tube distorted
# Isentropic Vortex
srun-serial ./laglos $options -of "${output_location}/iv" -m data/square5c0_vortex.mesh -p 5 -tf 2 -cfl 0.5 -rs 4 > "${output_files_dir}/out-iv" &

# IV distorted
# Noh
srun-serial ./laglos $options -of "${output_location}/noh" -m data/ref-square-c0.mesh -p 4 -tf 0.6 -cfl 0.25 -rs 5 > "${output_files_dir}/out-noh" &

# Noh distorted
# Sedov
srun-serial ./laglos $options -of "${output_location}/sedov" -m data/ref-square.mesh -p 6 -tf 0.9 -cfl 0.1 -rs 5 > "${output_files_dir}/out-sedov" &

# Saltzmann
srun-serial ./laglos $options -of "${output_location}/saltz" -m data/rectangle_saltzmann.mesh -p 7 -tf 0.6 -cfl 0.01 -rs 1 > "${output_files_dir}/out-saltz" &

# Triple Point
srun-serial ./laglos $options -of "${output_location}/tp" -m data/triple-point.mesh -p 12 -tf 5. -cfl 0.5 -rs 2 > "${output_files_dir}/out-tp" &


