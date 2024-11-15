"""
This python script file performs a convergence rates analsyis by running
a reading output files and composing a table from those values.

Example:
> python create_convergence_table.py "/Users/madisonsheridan/Workspace/Laglos/saved/convergence/temp_output/"
"""

import numpy as np
import re
import os
import sys
from tabulate import tabulate


# check command line inputs
assert len(sys.argv) == 2, "This file needs 1 input argument: directory, but " + str(len(sys.argv)) + " were given."
directory = str(sys.argv[1])
rate_precision = 6
# iterations = int(sys.argv[2])
# increment = int(sys.argv[3])
# plot_organization = False

# now we define the main function to be called at the end
def main():
    # comment out "run_simuations" if you only want to compute the errors
    vals = gather_vals()
    compute_rates(vals)


def gather_vals():
    vals = {'Processor_Runtime': [],
            'n_processes': [],
            'n_refinements': [],
            'n_Dofs': [],
            'h': [],
            # sv
            'sv_L1_Error': [],
            'sv_L1_Rates': [],
            'sv_L2_Error': [],
            'sv_L2_Rates': [],
            'sv_Linf_Error': [],
            'sv_Linf_Rates': [],
            # vel 
            'vel_L1_Error': [],
            'vel_L1_Rates': [],
            'vel_L2_Error': [],
            'vel_L2_Rates': [],
            'vel_Linf_Error': [],
            'vel_Linf_Rates': [],
            # ste
            'ste_L1_Error': [],
            'ste_L1_Rates': [],
            'ste_L2_Error': [],
            'ste_L2_Rates': [],
            'ste_Linf_Error': [],
            'ste_Linf_Rates': [],
            # composite
            'L1_Error': [],
            'L1_Rates': [],
            'L2_Error': [],
            'L2_Rates': [],
            'Linf_Error': [],
            'Linf_Rates': [],
            'mass_loss': [],
            'mass_loss_rates': [],
            'dt': [],
            'Endtime': []}
    for filename in sorted(os.listdir(directory)):
        f = os.path.join(directory, filename)
        with open(f) as fp:
            for cnt, ln in enumerate(fp):
                l = ln.strip().split()
                vals[l[0]].append(float(l[1]))

    return vals

def compute_rates(vals):

    ########## Use tabulate to create a formatted table for density rates ##########
    table = []
    for i in range(len(vals['h'])):
        if i == 0:
            table.append([vals['n_Dofs'][i], vals['sv_L1_Error'][i], "{---}", vals['sv_L2_Error'][i], "{---}",
                          vals['sv_Linf_Error'][i], "{---}"])
        else:
            if (vals['sv_L1_Error'][i-1] != 0):
                L1_rate = np.around(np.log(vals['sv_L1_Error'][i]/vals['sv_L1_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            else:
                L1_rate = "---"
            if (vals['sv_L2_Error'][i-1] != 0):
                L2_rate = np.around(np.log(vals['sv_L2_Error'][i]/vals['sv_L2_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            else:
                L2_rate = "---"
            if (vals['sv_Linf_Error'][i-1] != 0):
                Linf_rate = np.around(np.log(vals['sv_Linf_Error'][i]/vals['sv_Linf_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            else:
                Linf_rate = "---"
            # if vals['sv_L1_Error'][i-1] != 0:
            #    L1_rate = np.around(np.log(vals['sv_L1_Error'][i]/vals['sv_L1_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            # if vals['sv_L2_Error'][i-1] != 0:
            #    L2_rate = np.around(np.log(vals['sv_L2_Error'][i]/vals['sv_L2_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            # if vals['sv_Linf_Error'][i-1] != 0:
            #    Linf_rate = np.around(np.log(vals['sv_Linf_Error'][i]/vals['sv_Linf_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            # if vals['mass_loss'][i-1] > pow(10,-14):
            #     mass_loss_rate = np.around(np.log(vals['mass_loss'][i]/vals['mass_loss'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            table.append([vals['n_Dofs'][i], vals['sv_L1_Error'][i], L1_rate,
                     vals['sv_L2_Error'][i], L2_rate,
                     vals['sv_Linf_Error'][i], Linf_rate])
            # else:
            #     table.append([vals['n_Dofs'][i], vals['sv_L1_Error'][i], L1_rate,
            #               vals['sv_L2_Error'][i], L2_rate,
            #               vals['sv_Linf_Error'][i], Linf_rate,
            #               vals['mass_loss'][i], "{---}"])
            

    s_table = tabulate(table,
                       headers=["# dof", "sv L1 Error", "Rate", "sv L2 Error",
                                "Rate", "sv L-Inf Error", "Rate"],
                       tablefmt="latex")

    # Output table to console
    print("             ")
    print("Density Rates")
    print("             ")
    print(s_table)

    # Output table to txt file
    f = open("../saved/convergence/density_convergence_rates.txt", "w+")
    f.write(s_table)
    f.close()

    ########## Use tabulate to create a formatted table for velocity rates ##########
    table = []
    for i in range(len(vals['h'])):
        if i == 0:
            table.append([vals['n_Dofs'][i], vals['vel_L1_Error'][i], "{---}", vals['vel_L2_Error'][i], "{---}",
                          vals['vel_Linf_Error'][i], "{---}"])
        else:
            if (vals['vel_L1_Error'][i-1] != 0):
                L1_rate = np.around(np.log(vals['vel_L1_Error'][i]/vals['vel_L1_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            else:
                L1_rate = "---"
            if (vals['vel_L2_Error'][i-1] != 0):
                L2_rate = np.around(np.log(vals['vel_L2_Error'][i]/vals['vel_L2_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            else:
                L2_rate = "---"
            if (vals['vel_Linf_Error'][i-1] != 0):
                Linf_rate = np.around(np.log(vals['vel_Linf_Error'][i]/vals['vel_Linf_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            else:
                Linf_rate = "---"

            table.append([vals['n_Dofs'][i], vals['vel_L1_Error'][i], L1_rate,
                        vals['vel_L2_Error'][i], L2_rate,
                        vals['vel_Linf_Error'][i], Linf_rate])
            

    s_table = tabulate(table,
                       headers=["# dof", "vel L1 Error", "Rate", "vel L2 Error",
                                "Rate", "vel L-Inf Error", "Rate"],
                       tablefmt="latex")

    # Output table to console
    print("             ")
    print("Velocity Rates")
    print("             ")
    print(s_table)

    # Output table to txt file
    f = open("../saved/convergence/velocity_convergence_rates.txt", "w+")
    f.write(s_table)
    f.close()

    ########## Use tabulate to create a formatted table for specific total energy rates ##########
    table = []
    for i in range(len(vals['h'])):
        if i == 0:
            table.append([vals['n_Dofs'][i], vals['ste_L1_Error'][i], "{---}", vals['ste_L2_Error'][i], "{---}",
                          vals['ste_Linf_Error'][i], "{---}"])
        else:
            if (vals['ste_L1_Error'][i-1] != 0):
                L1_rate = np.around(np.log(vals['ste_L1_Error'][i]/vals['ste_L1_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            else:
                L1_rate = "---"
            if (vals['ste_L2_Error'][i-1] != 0):
                L2_rate = np.around(np.log(vals['ste_L2_Error'][i]/vals['ste_L2_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            else:
                L2_rate = "---"
            if (vals['ste_Linf_Error'][i-1] != 0):
                Linf_rate = np.around(np.log(vals['ste_Linf_Error'][i]/vals['ste_Linf_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            else:
                Linf_rate = "---"
            table.append([vals['n_Dofs'][i], vals['ste_L1_Error'][i], L1_rate,
                          vals['ste_L2_Error'][i], L2_rate,
                          vals['ste_Linf_Error'][i], Linf_rate])
            

    s_table = tabulate(table,
                       headers=["# dof", "ste L1 Error", "Rate", "ste L2 Error",
                                "Rate", "ste L-Inf Error", "Rate"],
                       tablefmt="latex")

    # Output table to console
    print("             ")
    print("Specific Total Energy Rates")
    print("             ")
    print(s_table)

    # Output table to txt file
    f = open("../saved/convergence/ste_convergence_rates.txt", "w+")
    f.write(s_table)
    f.close()

    ########## Use tabulate to create a formatted table for composite rates ##########
    table = []
    for i in range(len(vals['h'])):
        if i == 0:
            table.append([vals['n_Dofs'][i], vals['L1_Error'][i], "{---}", 
                        #   vals['L2_Error'][i], "{---}", vals['Linf_Error'][i], "{---}", 
                        #   vals['mass_loss'][i], "{---}"
                          ])
        else:
            if (vals['L1_Error'][i-1] != 0):
                L1_rate = np.around(np.log(vals['L1_Error'][i]/vals['L1_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            else:
                L1_rate = "---"
            if (vals['L2_Error'][i-1] != 0):
                L2_rate = np.around(np.log(vals['L2_Error'][i]/vals['L2_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            else:
                L2_rate = "---"
            if (vals['Linf_Error'][i-1] != 0):
                Linf_rate = np.around(np.log(vals['Linf_Error'][i]/vals['Linf_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            else:
                Linf_rate = "---"
            # L1_rate = np.around(np.log(vals['L1_Error'][i]/vals['L1_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            # L2_rate = np.around(np.log(vals['L2_Error'][i]/vals['L2_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            # Linf_rate = np.around(np.log(vals['Linf_Error'][i]/vals['Linf_Error'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
            if vals['mass_loss'][i-1] > pow(10,-14):
                mass_loss_rate = np.around(np.log(vals['mass_loss'][i]/vals['mass_loss'][i-1]) / np.log(vals['h'][i]/vals['h'][i-1]), decimals=rate_precision)
                table.append([vals['n_Dofs'][i], vals['L1_Error'][i], L1_rate,
                        #   vals['L2_Error'][i], L2_rate,
                        #   vals['Linf_Error'][i], Linf_rate,
                        #   vals['mass_loss'][i], mass_loss_rate
                          ])
            else:
                table.append([vals['n_Dofs'][i], vals['L1_Error'][i], L1_rate,
                        #   vals['L2_Error'][i], L2_rate,
                        #   vals['Linf_Error'][i], Linf_rate,
                        #   vals['mass_loss'][i], "{---}"
                          ])
            

    s_table = tabulate(table,
                       headers=["# dof", "L1 Error", "Rate", 
                              #   "L2 Error", "Rate", "L-Inf Error", "Rate", 
                              #   "Mass Loss", "Rate"
                                ],
                       tablefmt="latex")

    # Output table to console
    print("             ")
    print("Composite Rates")
    print("             ")
    print(s_table)

    # Output table to txt file
    f = open("../saved/convergence/composite_convergence_rates.txt", "w+")
    f.write(s_table)
    f.close()

# then we put main at the bottom to run everything
main()
