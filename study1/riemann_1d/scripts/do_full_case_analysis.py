# A script to do a complete analysis of an individual test case where multiple runs with varying
# dof and N have been done. It is assumed that for a given flux scheme, a case has been run with
# 200, 400 and 800 dofs in x-dir (nominally) and 12 dofs in y and z directions. It is also assumed
# that the N values are 1, 2, 3 and 5. The script generates the following plots
#
# 1. Visual comparison of results for given dof and varying N (outsourced)
# 2. Error vs dof for different values of N
# 3. Error vs wall time per time step for different values of N and dof
#
# All the required data will automatically be generated when 'do_full_individual_analysis.py' has
# been executed in a result directory.
#
# This script has to be run in the 'run' directory of a given case where individual simulation
# folders can be accessed. The generated plots are saved in '../plots' directory.
#
# Unlike 'do_full_individual_analysis.py', this is not designed to be called from an outer
# (shell/python) script and the data entries here have to be manually changed for every analysis.

import os
import subprocess

print("Doing full analysis in {}".format(os.getcwd()))
# directory where outsourced scripts lie
script_dir = "/home/vachan/Documents/Work/plens/study1/riemann_1d/scripts/"



N_values = [1,2,3,5]
dofs = [200, 400, 800]
actual_dof_values = [
    [100*2, 200*2, 400*2], # N=1
    [66*3, 134*3, 266*3], # N=2
    [50*4, 100*4, 200*4], # N=3
    [34*6, 66*6, 134*6], # N=5
]



# 1. Group plots showing results with different N for dixed dof
subprocess.run([
    "python3",
    script_dir + "plot_group.py",
    ".",
    str(200),
    "hllc",
    "result_logarithm",
    "comparison_data.csv",
    "Test 1, 200 dof, HLLC",
    "--save",
    "../plots",
    "dof200_12_12_hllc",
    "--size",
    "9",
    "6"
])
