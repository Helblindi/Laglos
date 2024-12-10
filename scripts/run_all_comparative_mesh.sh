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
options=" -mv 2 -ppd -greedy "
output_location="testing/20241210-mv2-ppd-greedy-em3/"
########

full_output_path="${results_dir}/${output_location}"
output_files_dir="${full_output_path}/output_files/"
mkdir "${results_dir}/${output_location}"
mkdir "${output_files_dir}"

echo $options

# Isentropic Vortex 
mkdir "${full_output_path}/iv-tf1"
mkdir "${full_output_path}/iv-tf10"
mkdir "${full_output_path}/iv-tf15"
srun-serial ./laglos $options -of "${output_location}/iv-tf1" -m data/square-vortex.mesh -p 5 -tf 1 -cfl 0.5 -rs 5 > "${output_files_dir}/out-ivtf1-rs5" &
srun-serial ./laglos $options -of "${output_location}/iv-tf1" -m data/square-vortex.mesh -p 5 -tf 1 -cfl 0.5 -rs 6 > "${output_files_dir}/out-ivtf1-rs6" &
srun-serial ./laglos $options -of "${output_location}/iv-tf10" -m data/square-vortex.mesh -p 5 -tf 10 -cfl 0.5 -rs 5 > "${output_files_dir}/out-ivtf10-rs5" &
srun-serial ./laglos $options -of "${output_location}/iv-tf10" -m data/square-vortex.mesh -p 5 -tf 10 -cfl 0.5 -rs 6 > "${output_files_dir}/out-ivtf10-rs6" &
srun-serial ./laglos $options -of "${output_location}/iv-tf15" -m data/square-vortex.mesh -p 5 -tf 15 -cfl 0.5 -rs 5 > "${output_files_dir}/out-ivtf15-rs5" &
srun-serial ./laglos $options -of "${output_location}/iv-tf15" -m data/square-vortex.mesh -p 5 -tf 15 -cfl 0.5 -rs 6 > "${output_files_dir}/out-ivtf15-rs6" &

# Sod tube
mkdir "${full_output_path}/sod-tube"
srun-serial ./laglos $options -of "${output_location}/sod-tube" -m data/shocktube.mesh -p 1 -tf .225 -cfl 0.5 -rs 3 > "${output_files_dir}/out-sod-tube-rs3" &
srun-serial ./laglos $options -of "${output_location}/sod-tube" -m data/shocktube.mesh -p 1 -tf .225 -cfl 0.5 -rs 4 > "${output_files_dir}/out-sod-tube-rs4" &

# Sod tube distorted
mkdir "${full_output_path}/sod-dis-tube"
srun-serial ./laglos $options -of "${output_location}/sod-dis-tube" -m data/distorted-tube-11.mesh -p 1 -tf .225 -cfl 0.5 -rs 3 > "${output_files_dir}/out-sod-dis-tube-rs3" &
srun-serial ./laglos $options -of "${output_location}/sod-dis-tube" -m data/distorted-tube-11.mesh -p 1 -tf .225 -cfl 0.5 -rs 4 > "${output_files_dir}/out-sod-dis-tube-rs4" &

# Sod radial
mkdir "${full_output_path}/sod-rad"
srun-serial ./laglos $options -of "${output_location}/sod-rad" -m data/ref-square-N15.mesh -p 13 -ti 0.01 -tf 0.2 -cfl 0.25 -rs 4 > "${output_files_dir}/out-sod-rad-rs4" &
srun-serial ./laglos $options -of "${output_location}/sod-rad" -m data/ref-square-N15.mesh -p 13 -ti 0.01 -tf 0.2 -cfl 0.25 -rs 5 > "${output_files_dir}/out-sod-rad-rs5" &

# Noh
mkdir "${full_output_path}/noh"
srun-serial ./laglos $options -of "${output_location}/noh" -m data/noh.mesh -p 4 -tf 0.6 -cfl 0.5 -rs 3 > "${output_files_dir}/out-noh-rs3" &
srun-serial ./laglos $options -of "${output_location}/noh" -m data/noh.mesh -p 4 -tf 0.6 -cfl 0.5 -rs 4 > "${output_files_dir}/out-noh-rs4" &

# Noh nonuniform
mkdir "${full_output_path}/noh-nonuniform"
srun-serial ./laglos $options -of "${output_location}/noh-nonuniform" -m data/noh-nonuniform.mesh -p 4 -tf 0.6 -cfl 0.5 -rs 1 > "${output_files_dir}/out-noh-rs1" &
srun-serial ./laglos $options -of "${output_location}/noh-nonuniform" -m data/noh-nonuniform.mesh -p 4 -tf 0.6 -cfl 0.5 -rs 2 > "${output_files_dir}/out-noh-rs2" &

# Sedov
mkdir "${full_output_path}/sedov"
srun-serial ./laglos $options -of "${output_location}/sedov" -m data/ref-square-N15.mesh -p 6 -tf 0.9 -cfl 1. -rs 4 > "${output_files_dir}/out-sedov-rs4" &
srun-serial ./laglos $options -of "${output_location}/sedov" -m data/ref-square-N15.mesh -p 6 -tf 0.9 -cfl 1. -rs 5 > "${output_files_dir}/out-sedov-rs5" &

# Saltzman
mkdir "${full_output_path}/saltzman"
srun-serial ./laglos $options -of "${output_location}/saltzman" -m data/rectangle_saltzmann.mesh -p 7 -tf 0.6 -cfl 0.01 -rs 1 > "${output_files_dir}/out-saltzman-rs1" &
srun-serial ./laglos $options -of "${output_location}/saltzman" -m data/rectangle_saltzmann.mesh -p 7 -tf 0.6 -cfl 0.01 -rs 2 > "${output_files_dir}/out-saltzman-rs2" &

# Triple Point
mkdir "${full_output_path}/tp"
srun-serial ./laglos $options -of "${output_location}/tp" -m data/triple-point.mesh -p 12 -tf 5. -cfl 0.5 -rs 2 > "${output_files_dir}/out-tp-rs2" &
srun-serial ./laglos $options -of "${output_location}/tp" -m data/triple-point.mesh -p 12 -tf 5. -cfl 0.5 -rs 3 > "${output_files_dir}/out-tp-rs3" &
srun-serial ./laglos $options -of "${output_location}/tp" -m data/triple-point.mesh -p 12 -tf 5. -cfl 0.5 -rs 4 > "${output_files_dir}/out-tp-rs4" &
srun-serial ./laglos $options -of "${output_location}/tp" -m data/triple-point.mesh -p 12 -tf 5. -cfl 0.5 -rs 5 > "${output_files_dir}/out-tp-rs5" &


