#!/bin/bash

file1="$1"
file2="$2"
num_tail_lines="${3:-11}"      # Lines to trim from the end
num_head_skip="${4:-10}"       # Lines to skip from the beginning

# Check if files exist
[ -f "$file1" ] || { echo "Missing file: $file1"; exit 1; }
[ -f "$file2" ] || { echo "Missing file: $file2"; exit 1; }

# Check if total lines are the same
if [ $(wc -l < "$file1") -ne $(wc -l < "$file2") ]; then
    echo "Files have different lengths. Test failed!"
    exit 1
fi

# Total lines minus head and tail
total_lines=$(wc -l < "$file1")
lines_to_keep=$((total_lines - num_head_skip - num_tail_lines))

if [ "$lines_to_keep" -le 0 ]; then
    echo "File is too short to compare after trimming."
    exit 1
fi

# Skip first N lines and trim last M lines
tail -n +"$((num_head_skip + 1))" "$file1" | head -n "$lines_to_keep" > file1_trimmed.txt
tail -n +"$((num_head_skip + 1))" "$file2" | head -n "$lines_to_keep" > file2_trimmed.txt

diff -u file1_trimmed.txt file2_trimmed.txt
exit $?
