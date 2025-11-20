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
mesh_file="${SOURCE_DIR}/data/elastic/human_leg_impact.mesh"
problem=53
final_time=.0012
cfl=0.5
max_steps=100
mv_option=2
refinements_serial=0

options="-m ${mesh_file} -p ${problem} "
options+="-tf ${final_time} -cfl ${cfl} "
options+="-mv ${mv_option} -ppd "
options+="-rs ${refinements_serial} "
options+="-ms ${max_steps} "
echo $options

# Run the executable with input file (adjust args as needed)
$EXECUTABLE $options > "${OUTPUT_FILE}"
