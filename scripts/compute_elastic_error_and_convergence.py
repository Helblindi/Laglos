import os
import glob
import numpy as np
import argparse
import re  # Import regex for extracting numbers
from tabulate import tabulate  # Import tabulate for LaTeX table generation

def extract_refinement_level(filename):
    """ Extracts refinement level from a filename like 'sv_06.csv' and computes step size """
    match = re.search(r'(\d+)', filename)  # Find the first number in the filename
    if match:
        level = int(match.group(1))  # Convert to integer
        return level, 2**-level  # Compute step size as h = 2^(-level)
    else:
        raise ValueError(f"Could not extract refinement level from {filename}")

def read_data(filename, column=0, has_header=False):
    """ 
    Reads a CSV or space-delimited TXT file and returns numerical values from the specified column.
    
    Parameters:
    - filename (str): Path to the file.
    - column (int): The column index to extract data from (0-based index).
    - has_header (bool): Whether the file has a header row to skip (True for approximation file).
    
    Returns:
    - np.array: Extracted numerical values.
    """
    delimiter = ',' if filename.endswith('.csv') else None  # None lets NumPy auto-detect spaces/tabs
    
    try:
        skip_rows = 1 if has_header else 0
        data = np.loadtxt(filename, delimiter=delimiter, usecols=[column], skiprows=skip_rows)
    except Exception as e:
        raise ValueError(f"Error reading file {filename}: {e}")
    
    return data

def interpolate_exact_solution(approx_x, exact_file, exact_col=0):
    """ Interpolates the exact solution onto the x-values of the approximation data. """
    exact_x = read_data(exact_file, 0, has_header=False)  # Read x-values from exact solution file
    exact_y = read_data(exact_file, exact_col, has_header=False)  # Read exact solution values
    
    interpolated_exact = np.interp(approx_x, exact_x, exact_y)  # Interpolate exact solution
    return interpolated_exact

def compute_L1_error(approx_file, exact_file, approx_col=0, exact_col=0):
    """ Computes the L1 error between approximation and exact solution. """
    approx_x = read_data(approx_file, 0, has_header=True)  # Read x-values from approx file
    approx = read_data(approx_file, approx_col, has_header=True)  # Approximation data (y-values)
    
    exact = interpolate_exact_solution(approx_x, exact_file, exact_col)
    
    if len(approx) != len(exact):
        raise ValueError("Mismatch in data size: Approximation and Exact Solution must have the same number of points.")
    
    L1_error = np.sum(np.abs(approx - exact)) / np.sum(np.abs(exact))
    return L1_error

def compute_convergence_order(errors, step_sizes):
    """ Computes the convergence order from L1 errors and step sizes. """
    orders = []
    for i in range(1, len(errors)):
        order = np.log(errors[i] / errors[i-1]) / np.log(step_sizes[i] / step_sizes[i-1])
        orders.append(order)
    return orders

def process_refinement_files(directory, exact_file, approx_col=0, exact_col=0):
    """ Processes the refinement files and computes the L1 error and convergence order. """
    approx_files = sorted(glob.glob(os.path.join(directory, "*.csv")))
    errors = []
    step_sizes = []
    
    for approx_file in approx_files:
        level, step_size = extract_refinement_level(os.path.basename(approx_file))
        L1_error = compute_L1_error(approx_file, exact_file, approx_col, exact_col)
        errors.append(L1_error)
        step_sizes.append(step_size)
    
    convergence_orders = compute_convergence_order(errors, step_sizes)
    
    return zip(approx_files, errors, step_sizes, convergence_orders)

def main():
    """ Main function to parse arguments and generate LaTeX table. """
    parser = argparse.ArgumentParser(description="Compute L1 error and convergence order for approximations.")
    parser.add_argument("directory", help="Directory containing the approximation files")
    parser.add_argument("exact_file", help="File containing the exact solution")
    parser.add_argument("approx_col", type=int, default=0, help="Column index for approximation data")
    parser.add_argument("exact_col", type=int, default=0, help="Column index for exact solution data")
    args = parser.parse_args()
    
    # Process the files and compute the results
    results = process_refinement_files(args.directory, args.exact_file, approx_col=args.approx_col, exact_col=args.exact_col)

    # Prepare data for tabulate
    table_data = []
    for file_name, error, step_size, order in results:
        step_size_str = f"{step_size:<10.3f}"
        error_str = f"{error:<2.3e}"
        order_str = f"{order:.3f}" if order != float('inf') else "N/A"
        
        # Add row to table_data
        table_data.append([f"{1. / step_size:<10.3f}", error_str, order_str])
    
    # Print LaTeX formatted table using tabulate
    headers = ['Num dofs', 'L1 Error', 'Convergence Order']
    latex_table = tabulate(table_data, headers=headers, tablefmt="latex", floatfmt=".3e")
    
    print(latex_table)

if __name__ == "__main__":
    main()
