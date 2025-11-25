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
mesh_file="${SOURCE_DIR}/data/shocktube.mesh"
problem=2
final_time=.225
cfl=0.5
mv_option=2
refinements_serial=0

options="-m ${mesh_file} -p ${problem} "
options+="-tf ${final_time} -cfl ${cfl} "
options+="-mv ${mv_option} -no-ppd "
options+="-rs ${refinements_serial} "
echo $options

# Run the executable with input file (adjust args as needed)
$EXECUTABLE $options > "${OUTPUT_FILE}"
