#! /bin/bash

# This is a script to run convergence analysis on Laglos using Release version
# These values are specific to Whistler
laglos_dir="/Users/madisonsheridan/Workspace/Laglos"
scripts_dir="${laglos_dir}/scripts"
bin_dir="${laglos_dir}/build"
create_convergence="${scripts_dir}/create_convergence_table_sv.py"
results_dir="${bin_dir}/results/tests"


cd ${bin_dir}

####### Global test params
mv_option=2
fv_option=2
base_options="-mv ${mv_option} -fv ${fv_option} -ppd "

###################################
####### DISTORTED SOD TUBE TEST
###################################
mesh_file="data/distorted-tube-11.mesh"
problem=1
final_time=.225
cfl=0.5
mv_iter_n=0
output_location="tests/sod-dis"
########

options+="${base_options} -m ${mesh_file} -p ${problem} "
options+="-tf ${final_time} -cfl ${cfl} "
options+="-of ${output_location} "

./Laglos -rs 0 $options
./laglos -rs 1 $options

###################################
####### NOH TEST
###################################
# mesh_file="data/noh.mesh"
# problem=4
# final_time=0.6
# cfl=0.5
# output_location="tests/noh"
# ########

# options="${base_options} -m ${mesh_file} -p ${problem} "
# options+="-tf ${final_time} -cfl ${cfl} "
# options+="-of ${output_location} "

# ./laglos -rs 0 $options
# ./laglos -rs 1 $options

# ###################################
# ####### ISENTROPIC VORTEX TEST
# ###################################
# mesh_file="data/square-vortex.mesh"
# problem=5
# final_time=1
# cfl=0.5
# output_location="tests/iv-tf1"
# ########

# options="${base_options} -m ${mesh_file} -p ${problem} "
# options+="-tf ${final_time} -cfl ${cfl} "
# options+="-of ${output_location} "

# ./laglos -rs 1 $options
# ./laglos -rs 2 $options
# ./laglos -rs 3 $options

# ###################################
# ####### SEDOV TEST
# ###################################
# mesh_file="data/ref-square-N15.mesh"
# problem=6
# final_time=0.9
# cfl=1.
# output_location="tests/sedov"
# ########

# options="${base_options} -m ${mesh_file} -p ${problem} "
# options+="-tf ${final_time} -cfl ${cfl} "
# options+="-of ${output_location} "

# ./laglos -rs 0 $options
# ./laglos -rs 1 $options

# # Create convergence tables for tests
python3 ${create_convergence} ${results_dir}/sod-dis/convergence/temp_output
# python3 ${create_convergence} ${results_dir}/noh/convergence/temp_output
# python3 ${create_convergence} ${results_dir}/iv-tf1/convergence/temp_output
# python3 ${create_convergence} ${results_dir}/sedov/convergence/temp_output
