#!/bin/bash

file1="$1"
file2="$2"
num_lines="${3:-11}"  # Default to ignoring last 11 lines

# Check if files exist
[ -f "$file1" ] || { echo "Missing file: $file1"; exit 1; }
[ -f "$file2" ] || { echo "Missing file: $file2"; exit 1; }

# Check if total lines are the same
if [ $(wc -l < "$file1") -ne $(wc -l < "$file2") ]; then
    echo "Files have different lengths. Test failed!"
    exit 1
fi


# Count total lines
total_lines=$(wc -l < "$file1")
lines_to_keep=$((total_lines - num_lines))

# If fewer lines than lines_to_ignore, skip comparison
if [ "$lines_to_keep" -le 0 ]; then
    echo "File is too short to compare after trimming."
    exit 1
fi

# Trim and compare
sed "${lines_to_keep}q" "$file1" > file1_trimmed.txt
sed "${lines_to_keep}q" "$file2" > file2_trimmed.txt

diff -u file1_trimmed.txt file2_trimmed.txt
exit $?
