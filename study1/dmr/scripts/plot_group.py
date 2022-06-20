# A script to do line plots of extracted data. This script can be executed from any location as
# long as the `run_dir` and `res_dirs` variables are set correctly.

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 14,
    "axes.formatter.limits": [-2,2]
})
import pandas as pd

# the common directory for all runs
run_dir = "/home/vachan/Documents/Work/plens/study1/dmr/run"
N_values = [1,2,3,5]
flux = "chandrashekhar" # flux scheme used, generally this is a suffix to the case dir
res_dir_name = "result" # name of result dir within case dir
line_loc_suffixes = ["x2.5", "y0.3"] # the suffixes for the different number of lines sampled
plot_dir = "/home/vachan/Documents/Work/plens/study1/dmr/plots" # loc for saving plots

# set the data file names
data_filenames = pd.DataFrame(data="", index=N_values, columns=line_loc_suffixes)
for N in N_values:
    res_dir = "{}/res_60_N{}_{}/{}/".format(run_dir, N, flux, res_dir_name)
    for s in line_loc_suffixes:
        data_filenames.loc[N, s] = res_dir+ "line_data_{}.csv".format(s)

# set the axes variables for different lines
xaxis_vars = pd.Series(data=["Points1", "Points0"], index=line_loc_suffixes)
xaxis_labels = pd.Series(data=[r"$y$", r"$x$"], index=line_loc_suffixes)
yaxis_vars = pd.Series(data=["rhov", "T"], index=line_loc_suffixes)
yaxis_labels = pd.Series(data=[r"$\rho v$", r"$T$"], index=line_loc_suffixes)

N_colors = pd.Series(data=["blue", "green", "red", "magenta"], index=N_values)
N_linestyles = pd.Series(
    data=[
        (0, ()), # solid line
        (0, (1,1)), # dotted
        (0, (5,1)), # dense dash
        (0, (5,2,1,2)) # dash dot
    ],
    index=N_values
)
for line_suffix in line_loc_suffixes:
    fig, ax = plt.subplots(1,1)
    for N in N_values:
        data = np.genfromtxt(data_filenames.loc[N, line_suffix], delimiter=",", names=True)
        ax.plot(
            data[xaxis_vars.loc[line_suffix]],
            data[yaxis_vars.loc[line_suffix]],
            label=r"$N={}$".format(N),
            ls=N_linestyles.loc[N],
            marker="",
            c=N_colors.loc[N],
            alpha=0.75
        )
    ax.set_xlabel(xaxis_labels.loc[line_suffix])
    ax.set_ylabel(yaxis_labels.loc[line_suffix])
    ax.grid()
    ax.legend()
    fig.tight_layout(pad=0.25)
    plt.show()
    for fmt in ["png", "pdf"]:
        fig_filename = "{}/line_{}.{}".format(plot_dir, line_suffix, fmt)
        fig.savefig(fig_filename, format=fmt)
        print("Saved file {}".format(fig_filename))
