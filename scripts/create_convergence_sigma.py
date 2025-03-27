import numpy as np
import os
import sys
from tabulate import tabulate

# check command line inputs
assert len(sys.argv) == 2, "This file needs 1 input argument: directory, but " + str(len(sys.argv)) + " were given."
directory = str(sys.argv[1])
rate_precision = 6

def main():
    vals = gather_vals()
    compute_rates(vals)

def gather_vals():
    vals = {
        'Processor_Runtime': [], 'n_processes': [], 'n_refinements': [], 'n_Dofs': [], 'h': [],
        'sv_L1_Error': [], 'sv_L2_Error': [], 'sv_Linf_Error': [],
        'vel_L1_Error': [], 'vel_L2_Error': [], 'vel_Linf_Error': [],
        'ste_L1_Error': [], 'ste_L2_Error': [], 'ste_Linf_Error': [],
        'L1_Error': [], 'L2_Error': [], 'Linf_Error': [], 'mass_loss': [],
        'mass_loss_rates': [], 'dt': [], 'Endtime': []
    }

    sigma_keys = ['sigma_L1_Error', 'sigma_L2_Error', 'sigma_Linf_Error']  # NEW
    sigma_exists = {key: False for key in sigma_keys}  # NEW

    for filename in sorted(os.listdir(directory)):
        f = os.path.join(directory, filename)
        with open(f) as fp:
            for cnt, ln in enumerate(fp):
                l = ln.strip().split()
                key = l[0]
                if key in vals:
                    vals[key].append(float(l[1]))
                elif key in sigma_keys:  # NEW: Only add if sigma errors exist
                    if key not in vals:
                        vals[key] = []
                    vals[key].append(float(l[1]))
                    sigma_exists[key] = True  # NEW: Track presence of sigma errors

    vals['sigma_exists'] = any(sigma_exists.values())  # NEW: Flag if sigma exists
    return vals

def compute_rates(vals):
    def compute_error_rates(error_list, h_list):
        rates = ["{---}"]
        for i in range(1, len(error_list)):
            if error_list[i - 1] != 0:
                rate = np.around(np.log(error_list[i] / error_list[i - 1]) / np.log(h_list[i] / h_list[i - 1]), decimals=rate_precision)
                rates.append(rate)
            else:
                rates.append("{---}")
        return rates

    def print_table(errors, label):
        table = []
        h_list = vals['h']
        for i in range(len(h_list)):
            row = [vals['n_Dofs'][i]]
            for err_type in errors:
                row.append(vals[err_type][i])
                if i == 0:
                    row.append("{---}")
                else:
                    row.append(compute_error_rates(vals[err_type], h_list)[i])
            table.append(row)

        headers = ["# dof"] + [f"{e.split('_')[0]} {e.split('_')[1]}" for e in errors for _ in (0, 1)]
        s_table = tabulate(table, headers=headers, tablefmt="latex")

        print(f"\n{label} Rates\n")
        print(s_table)

        with open(f"../saved/convergence/{label.lower().replace(' ', '_')}_convergence_rates.txt", "w+") as f:
            f.write(s_table)

    # Print tables for density, velocity, and specific total energy
    print_table(['sv_L1_Error', 'sv_L2_Error', 'sv_Linf_Error'], "Density")
    print_table(['vel_L1_Error', 'vel_L2_Error', 'vel_Linf_Error'], "Velocity")
    print_table(['ste_L1_Error', 'ste_L2_Error', 'ste_Linf_Error'], "Specific Total Energy")

    # NEW: Print sigma table **only if sigma data exists**
    if vals.get('sigma_exists', False):
        print_table(['sigma_L1_Error', 'sigma_L2_Error', 'sigma_Linf_Error'], "Sigma")

# Run everything
main()
