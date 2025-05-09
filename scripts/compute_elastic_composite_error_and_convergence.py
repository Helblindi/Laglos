import os
import glob
import numpy as np
import pandas as pd
import argparse
import re
from tabulate import tabulate


def file_has_header(filepath):
    """Determines if the file has a header by attempting to parse the first line as floats."""
    with open(filepath, 'r') as f:
        first_line = f.readline().strip().split(',')
        try:
            [float(item) for item in first_line]
            return False  # All values numeric → no header
        except ValueError:
            return True   # Non-numeric value found → header present


def extract_refinement_level(filename):
    """Extracts refinement level from a filename and computes the step size as 2^-level."""
    match = re.search(r'(\d+)', filename)
    if match:
        level = int(match.group(1))
        return level, 2**-level
    else:
        raise ValueError(f"Could not extract refinement level from {filename}")


def read_data(filename, column=0, has_header=False):
    """Reads a CSV file and extracts a single column of numeric data."""
    try:
        df = pd.read_csv(filename, header=0 if has_header else None)
        return df.iloc[:, column].astype(float).to_numpy()
    except Exception as e:
        raise ValueError(f"Error reading file {filename}: {e}")


def interpolate_exact_solution(approx_x, exact_file, exact_col=0):
    """Interpolates the exact solution onto the x-values of the approximation data."""
    _has_header = file_has_header(exact_file)
    exact_x = read_data(exact_file, 0, has_header=_has_header)
    exact_y = read_data(exact_file, exact_col, has_header=_has_header)
    return np.interp(approx_x, exact_x, exact_y)


def compute_L1_error(approx_file, exact_file, approx_cols, exact_cols):
    """Computes the composite L1 error between approximation and exact solution across variables."""
    approx_x = read_data(approx_file, 0, has_header=True)
    approx = [read_data(approx_file, col, has_header=True) for col in approx_cols]
    exact = [interpolate_exact_solution(approx_x, exact_file, col) for col in exact_cols]
    return sum(np.sum(np.abs(a - e)) / np.sum(np.abs(e)) for a, e in zip(approx, exact))


def compute_L2_error(approx_file, exact_file, approx_cols, exact_cols):
    """Computes the composite L2 error between approximation and exact solution across variables."""
    approx_x = read_data(approx_file, 0, has_header=True)
    approx = [read_data(approx_file, col, has_header=True) for col in approx_cols]
    exact = [interpolate_exact_solution(approx_x, exact_file, col) for col in exact_cols]
    return sum(np.sqrt(np.sum((a - e) ** 2)) / np.sqrt(np.sum(e ** 2)) for a, e in zip(approx, exact))


def compute_convergence_order(errors, step_sizes):
    """Computes convergence order using the formula log(Ei/Ei-1) / log(hi/hi-1)."""
    return [np.log(errors[i] / errors[i - 1]) / np.log(step_sizes[i] / step_sizes[i - 1]) for i in range(1, len(errors))]


def process_refinement_files(directory, exact_file, approx_cols, exact_cols):
    """Processes refinement files to compute L1 errors and convergence orders."""
    approx_files = sorted(glob.glob(os.path.join(directory, "*.csv")))
    errors = []
    step_sizes = []
    for approx_file in approx_files:
        level, step_size = extract_refinement_level(os.path.basename(approx_file))
        errors.append(compute_L1_error(approx_file, exact_file, approx_cols, exact_cols))
        step_sizes.append(step_size)
    convergence_orders = compute_convergence_order(errors, step_sizes)
    return errors, step_sizes, convergence_orders


def main():
    """Main function to parse arguments, compute errors, and print a LaTeX-formatted table."""
    parser = argparse.ArgumentParser(description="Compute composite L1 error and convergence order for approximations.")
    parser.add_argument("directory", help="Directory containing the approximation files")
    parser.add_argument("exact_file", help="File containing the exact solution")
    parser.add_argument("--approx_cols", type=int, nargs='+', default=[2, 3, 5])
    parser.add_argument("--exact_cols", type=int, nargs='+', default=[13, 18, 15])
    args = parser.parse_args()

    errors, step_sizes, convergence_orders = process_refinement_files(args.directory, args.exact_file, args.approx_cols, args.exact_cols)

    table_data = []
    for i, error in enumerate(errors):
        order_str = "---" if i == 0 else f"{convergence_orders[i - 1]:<0.2f}"
        num_dofs = int(1 / step_sizes[i])
        table_data.append([num_dofs, error, order_str])

    headers = ['Num dofs', 'Composite L1 Error', 'Convergence Order']
    print(tabulate(table_data, headers=headers, tablefmt="latex", floatfmt=(".0f", "1.3e")))


if __name__ == "__main__":
    main()
