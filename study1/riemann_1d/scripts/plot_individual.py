# Plots individual results by comparing with analytical solution.

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14
plt.rcParams["mathtext.fontset"] = "dejavuserif"
from ic_data import *
from shockTubeSoln import shockTubeSoln

parser = argparse.ArgumentParser(
    description = "A script to compare individual results with reference solution. Also prints "
        + "some error metrics. The error is calculated evaluating the difference at the "
        + "simulation data points provided in 'sim_data_filename'."
)
parser.add_argument(
    "sim_data_filename",
    help="The file containing extracted simulation data for comparison with reference solution. "
        + "This file can be generated using 'extract_data.py'. Full path has to be provided here, "
        +"result directory will be extracted from this argument and used for saving plots. Assumed "
        +"to be in csv format with 1st row headers."
)
parser.add_argument(
    "test",
    help="The test case being setup. Options are: {}.".format(test_options)
)
parser.add_argument(
    "dia_loc",
    help="The exact diaphragm location in the simulation. In many cases, since simulations set "
        + "diaphragm "
        + "at cell interface, it may not be possible to exactly set the diaphragm at required "
        + "location.",
    type=float
)
parser.add_argument(
    "plot_filename",
    help="Base file name of the plot to be saved. The plot is saved in '.png' and '.pdf' formats "
        + "in the result directory extracted from 'sim_data_filename'."
)
args = parser.parse_args()

print("Reading simulation data from {}".format(args.sim_data_filename))
sim_data = np.genfromtxt(args.sim_data_filename, delimiter=",", names=True)
ex_data = shockTubeSoln(
    p_left[args.test], p_right[args.test],
    p_left[args.test]/(287*rho_left[args.test]), p_right[args.test]/(287*rho_right[args.test]),
    u_left[args.test], u_right[args.test],
    sim_data["Points0"]-args.dia_loc,
    end_times[args.test]
)

fig, ax = plt.subplots(1,1)
ax.plot(sim_data["Points0"], ex_data[:,4], "b-", label="Exact")
ax.plot(sim_data["Points0"], sim_data["rho"], "r-", alpha=0.75, label="Simulation")
ax.legend(loc="best")
ax.grid()
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$\rho$")
fig.tight_layout()

plt.show()

res_dir = os.path.dirname(args.sim_data_filename) + "/"
for fmt in ["png", "pdf"]:
    full_plot_filename = res_dir + args.plot_filename + "." + fmt
    fig.savefig(full_plot_filename, format=fmt)
    print("Written file {}".format(full_plot_filename))
full_data_filename = res_dir + args.plot_filename + "_data.csv"
np.savetxt(
    full_data_filename,
    np.transpose(np.vstack((sim_data["Points0"], ex_data[:,4], sim_data["rho"]))),
    delimiter=","
)
print("Written file {}".format(full_data_filename))



error_vec = sim_data["rho"] - ex_data[:,4]
errors_content = """Error norm order, Error
1, {}
2, {}
inf, {}
""".format(
    np.linalg.norm(error_vec, ord=1)/np.linalg.norm(ex_data[:,4], ord=1),
    np.linalg.norm(error_vec, ord=2)/np.linalg.norm(ex_data[:,4], ord=2),
    np.linalg.norm(error_vec, ord=np.inf)/np.linalg.norm(ex_data[:,4], ord=np.inf)
)
print("Errors:")
print(errors_content)
full_errors_filename = res_dir + "errors.csv"
errors_file = open(full_errors_filename, "w")
errors_file.write(errors_content)
errors_file.close()
print("Written to file {}".format(full_errors_filename))
