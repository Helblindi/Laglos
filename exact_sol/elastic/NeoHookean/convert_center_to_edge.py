import sys

def add_offset(input_file, output_file, offset=0.0005):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            values = line.split()
            new_values = [str(float(value) + offset) for value in values]
            outfile.write(" ".join(new_values) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 add_offset.py <input_file> <output_file>")
        sys.exit(1)
    
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    
    add_offset(input_filename, output_filename)