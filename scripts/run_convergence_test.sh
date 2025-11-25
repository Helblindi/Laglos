#! /bin/bash
# Resolve source and build dirs (prefer CMake-configured values if script was configured with configure_file)
LAGLOS_SOURCE_DIR="@CMAKE_SOURCE_DIR@"
LAGLOS_BINARY_DIR="@CMAKE_BINARY_DIR@"

# If not configured by CMake (placeholders remained), derive from script location
if [ "@CMAKE_SOURCE_DIR@" = "@""CMAKE_SOURCE_DIR""@" ]; then
   script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
   LAGLOS_SOURCE_DIR="$(cd "${script_dir}/.." && pwd)"
   if [ -d "${LAGLOS_SOURCE_DIR}/build" ]; then
      LAGLOS_BINARY_DIR="${LAGLOS_SOURCE_DIR}/build"
   else
      LAGLOS_BINARY_DIR="${LAGLOS_SOURCE_DIR}"
   fi
fi

export LAGLOS_SOURCE_DIR LAGLOS_BINARY_DIR
scripts_dir="${LAGLOS_SOURCE_DIR}/scripts"

# Fix OpenMPI shared memory issues on macOS
# Use TCP instead of shared memory for single-node communication
export OMPI_MCA_btl=^sm

echo "Laglos source dir: ${LAGLOS_SOURCE_DIR}"
echo "Laglos binary dir: ${LAGLOS_BINARY_DIR}"

cd $LAGLOS_BINARY_DIR


# EDIT THESE PARAMS
# results_dir="/home/sheridanm/scratch/Laglos-results"
results_dir="${LAGLOS_BINARY_DIR}/results"
mesh_file="../data/ref-segment.mesh"
problem=2
final_time=.225
cfl=0.5
mv_option=2
fv_option=0
# mv_it_op=2
# mv_iter_n=2
# mm_visc=0.
output_location="testing/sod"
output_file_prefix="${results_dir}/${output_location}/out-sod-r"
elastic_op=0
MIN_RS=3
MAX_RS=8
########

# Create output directory
mkdir -p "${results_dir}/${output_location}"

# Set options
options="-m ${mesh_file} -p ${problem} "
options+="-tf ${final_time} -cfl ${cfl} "
options+="-mv ${mv_option} -fv ${fv_option} "
# options+="-mv-it-op ${mv_it_op} -mv-iter-n ${mv_iter_n} " 
# options+="-mmv ${mm_visc} "
options+="-of ${results_dir}/${output_location} "
options+="-pview -vs 5 "
options+="-ue ${elastic_op} "

echo "Options provided: $options"

# Run with varying rs
min_rs="${MIN_RS:-1}"
max_rs="${MAX_RS:-8}"
for rs in $(seq "$min_rs" "$max_rs"); do
   ## Command to run in serial on given machine
   # Whistler example:
   # srun-serial -n 1 ./Laglos -rs "$rs" $options > "${output_file_prefix}${rs}" &
   # Local example:
   ./Laglos -rs "$rs" $options > "${output_file_prefix}${rs}" &
done

# Lastly, copy this script to the output_file location
# for documentation on run parameters
cp "${scripts_dir}/run_convergence_test.sh" "${results_dir}/${output_location}/." 
