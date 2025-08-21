#! /bin/bash

SOURCE_DIR="$1"
OUTPUT_FILE="$2"

# Full path to the executable
EXECUTABLE="${SOURCE_DIR}/build/Laglos"

# Check that the executable exists
if [[ ! -x "$EXECUTABLE" ]]; then
    echo "Executable not found or not executable: $EXECUTABLE"
    exit 1
fi

# test parameters
mesh_file="${SOURCE_DIR}/data/ref-square.mesh"
problem=0
final_time=.1
cfl=0.5
mv_option=2
refinements_serial=3
order_t=1
order_k=2
visc_option=3
solver=13

options="-m ${mesh_file} -p ${problem} "
options+="-tf ${final_time} -cfl ${cfl} "
options+="-mv ${mv_option} "
options+="-rs ${refinements_serial} -ppd "
options+="-ot ${order_t} -ok ${order_k} "
options+="-visc ${visc_option} -s ${solver} "
echo $options

# Run the executable with input file (adjust args as needed)
$EXECUTABLE $options > "${OUTPUT_FILE}"
