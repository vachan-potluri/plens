# Plots individual results by comparing with reference solution.

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14
plt.rcParams["mathtext.fontset"] = "dejavuserif"

parser = argparse.ArgumentParser(
    description = "A script to compare individual results with reference solution."
)
parser.add_argument(
    "sim_data_filename",
    help="The file containing extracted simulation data for comparison with reference solution. "
        + "This file can be generated using 'extract_data.py'. Full path has to be provided here, "
        +"result directory will be extracted from this argument and used for saving plots."
)
parser.add_argument(
    "ref_data_filename",
    help="The file containing reference data for comparison with simulation. "
        + "Full path is required."
)
parser.add_argument(
    "plot_filename",
    help="File name of the plot to be saved. The plot is saved in '.png' and '.pdf' formats in "
        + "the result directory extracted from 'sim_data_filename'."
)
args = parser.parse_args()

print("Reading simulation data from {}".format(args.sim_data_filename))
sim_data = np.genfromtxt(args.sim_data_filename, delimiter=",", names=True)
print("Reading reference data from {}".format(args.ref_data_filename))
ref_data = np.genfromtxt(args.ref_data_filename, delimiter=",")

fig, ax = plt.subplots(1,1)
ax.plot(ref_data[:,0], ref_data[:,1], "b-", label="Reference\n(Shu & Osher, 1989)")
ax.plot(sim_data["Points0"], sim_data["rho"], "ro", markersize=2, label="Simulation")
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
