# Plots individual results by comparing with reference solution.

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 14
plt.rcParams["mathtext.fontset"] = "dejavuserif"
from scipy.interpolate import interp1d

parser = argparse.ArgumentParser(
    description = "A script to compare individual results with reference solution. Also prints "
        + "some error metrics. The error is calculated evaluating the difference at the "
        + "simulation data points provided in 'sim_data_filename'. Reference solution is "
        + "interpolated in piecewise cubic manner for this purpose."
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
    help=" Base file name of the plot to be saved. The plot is saved in '.png' and '.pdf' formats "
        + "in the result directory extracted from 'sim_data_filename'."
)
args = parser.parse_args()

print("Reading simulation data from {}".format(args.sim_data_filename))
sim_data = np.genfromtxt(args.sim_data_filename, delimiter=",", names=True)
print("Reading reference data from {}".format(args.ref_data_filename))
ref_data = np.genfromtxt(args.ref_data_filename, delimiter=",")

fig, ax = plt.subplots(1,1)
ax.plot(ref_data[:,0], ref_data[:,1], "b-", label="Reference\n(Shu & Osher, 1989)")
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

x_ref = ref_data[:,0]
rho_ref = ref_data[:,1]
ref_func = interp1d(x_ref, rho_ref, kind="cubic")
x_sim = sim_data["Points0"]
mask = (x_sim > np.min(x_ref)) & (x_sim < np.max(x_ref))
rho_ref_evaluated = ref_func(x_sim[mask])
error_vec = rho_ref_evaluated - sim_data["rho"][mask]
errors_content = """Error norm order, Error
1, {}
2, {}
inf, {}
""".format(
    np.linalg.norm(error_vec, ord=1)/np.linalg.norm(rho_ref_evaluated, ord=1),
    np.linalg.norm(error_vec, ord=2)/np.linalg.norm(rho_ref_evaluated, ord=2),
    np.linalg.norm(error_vec, ord=np.inf)/np.linalg.norm(rho_ref_evaluated, ord=np.inf)
)
print("Errors:")
print(errors_content)
full_errors_filename = res_dir + "errors.csv"
errors_file = open(full_errors_filename, "w")
errors_file.write(errors_content)
errors_file.close()
print("Written to file {}".format(full_errors_filename))
