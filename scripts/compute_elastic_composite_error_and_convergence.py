import os
import glob
import numpy as np
import argparse
import re
from tabulate import tabulate

def extract_refinement_level(filename):
    """ Extracts refinement level from a filename like 'sv_06.csv' and computes step size """
    match = re.search(r'(\d+)', filename)  # Find the first number in the filename
    if match:
        level = int(match.group(1))  # Convert to integer
        return level, 2**-level  # Compute step size as h = 2^(-level)
    else:
        raise ValueError(f"Could not extract refinement level from {filename}")

def read_data(filename, column=0, has_header=False):
    """ Reads a CSV or space-delimited TXT file and returns numerical values from the specified column. """
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

def compute_L1_error(approx_file, exact_file, approx_cols, exact_cols):
    """ Computes the composite L1 error between approximation and exact solution across multiple variables. """
    approx_x = read_data(approx_file, 0, has_header=True)  # Read x-values from approx file
    approx = [read_data(approx_file, col, has_header=True) for col in approx_cols]  # Approximation data (y-values)
    
    exact = [interpolate_exact_solution(approx_x, exact_file, col) for col in exact_cols]  # Exact solution for all columns
    
    # Compute the composite L1 error (sum of relative errors across all variables)
    _errors = (np.sum(np.abs(a - e)) / np.sum(np.abs(e)) for a, e in zip(approx, exact))
    composite_error = sum(np.sum(np.abs(a - e)) / np.sum(np.abs(e)) for a, e in zip(approx, exact))
    return composite_error

def compute_convergence_order(errors, step_sizes):
    """ Computes the convergence order using the formula: log(Ei / Ei-1) / log(hi / hi-1) """
    orders = []
    for i in range(1, len(errors)):
        E_i = errors[i]
        E_i_1 = errors[i-1]
        h_i = step_sizes[i]
        h_i_1 = step_sizes[i-1]
        order = np.log(E_i / E_i_1) / np.log(h_i / h_i_1)
        orders.append(order)
    return orders

def process_refinement_files(directory, exact_file, approx_cols, exact_cols):
    """ Processes the refinement files and computes the composite L1 error and convergence order. """
    approx_files = sorted(glob.glob(os.path.join(directory, "*.csv")))
    errors = []
    step_sizes = []
    
    for approx_file in approx_files:
        level, step_size = extract_refinement_level(os.path.basename(approx_file))
        L1_error = compute_L1_error(approx_file, exact_file, approx_cols, exact_cols)
        errors.append(L1_error)
        step_sizes.append(step_size)
    
    convergence_orders = compute_convergence_order(errors, step_sizes)

    return errors, step_sizes, convergence_orders

def main():
    """ Main function to parse arguments and generate LaTeX table. """
    parser = argparse.ArgumentParser(description="Compute composite L1 error and convergence order for approximations.")
    parser.add_argument("directory", help="Directory containing the approximation files")
    parser.add_argument("exact_file", help="File containing the exact solution")
    parser.add_argument("--approx_cols", type=int, nargs='+', default=[0, 1, 2], help="Column indices for approximation data (density, velocity, energy)")
    parser.add_argument("--exact_cols", type=int, nargs='+', default=[1, 2, 3], help="Column indices for exact solution data (density, velocity, energy)")
    args = parser.parse_args()
    
    # Process the files and compute the results
    errors, step_sizes, convergence_orders = process_refinement_files(args.directory, args.exact_file, approx_cols=args.approx_cols, exact_cols=args.exact_cols)

    # Prepare data for tabulate
    table_data = []
    
    for i, error in enumerate(errors):

        error_str = f"{error:<1.3e}"

        # For the first row, do not compute convergence order, display '{---}'
        if i == 0:
            order_str = "---"
        else:
            order = convergence_orders[i-1]
            order_str = f"{order:<0.2f}"
        
        # Add row to table data
        num_dofs = int(1 / step_sizes[i])
        table_data.append([num_dofs, error, order_str])
    
    # Print LaTeX formatted table using tabulate
    headers = ['Num dofs', 'Composite L1 Error', 'Convergence Order']
    latex_table = tabulate(table_data, headers=headers, tablefmt="latex", floatfmt=(".0f", "1.3e"))
    
    print(latex_table)

if __name__ == "__main__":
    main()
