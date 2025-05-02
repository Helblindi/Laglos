import sys

def split_columns(input_file, columns):
    with open(input_file, 'r') as f:
        lines = [line.strip().split() for line in f]

    # Validate the maximum column index
    max_index = max(columns)
    for i, line in enumerate(lines):
        if len(line) <= max_index:
            print(f"Warning: Line {i+1} has only {len(line)} columns, but index {max_index} was requested.")
            continue  # Skip this line to avoid IndexError

        for col in columns:
            with open(f'column_{col}.txt', 'a') as out_f:
                out_f.write(line[col] + '\n')

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python col_parser.py <input_file> <col1> <col2> ...")
        sys.exit(1)

    input_file = sys.argv[1]
    columns = [int(c) for c in sys.argv[2:]]  # Convert column indices to integers

    split_columns(input_file, columns)
