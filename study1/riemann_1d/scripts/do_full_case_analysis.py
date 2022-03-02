# A script to do a complete analysis of an individual test case where multiple runs with varying
# dof and N have been done. It is assumed that for a given flux scheme, a case has been run with
# 200, 400 and 800 dofs in x-dir (nominally) and 12 dofs in y and z directions. It is also assumed
# that the N values are 1, 2, 3 and 5. The script generates the following plots
#
# 1. Visual comparison of results for given dof and varying N (outsourced)
# 2. Error vs dof for different values of N
# 3. Error vs cpu time per time step for different values of N and across dofs
# 4. Visual comparison of results for all dofs and all values of N (outsourced)
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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter, MultipleLocator
from scipy.stats import linregress
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams["font.size"] = 10 # ineffective
plt.rcParams["axes.formatter.limits"] = [-2,2]
plt.rcParams["axes.formatter.use_mathtext"] = True



def find_nearest(array, value):
    # Finds and returns an `array` element closest to `value`
    # Source: https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]



def mysavefig(fig, save_dir, base_name):
    # Saves the figure handle
    for fmt in ["pdf", "png"]:
        filename = "{}/{}.{}".format(save_dir, base_name, fmt)
        fig.savefig(filename, format=fmt)
        print("Saved file {}".format(filename))



def plot_convergence_rate(ax, x, y, loc="up"):
    # From the given `x`-`y` line data, this function calculates nearest convergence rate in
    # interval of 0.5 and plots it on the `ax` provided. The `loc` parameter is used to determine
    # whether the rate will be shown above or below the line data (this function doesn't plot the
    # line plot itself)
    # For `ax` with multiple line plots, the bottom-most or top-most line data can be passed to
    # this function with `loc="below"` or `loc="up"`
    x = np.array(x)
    y = np.array(y)
    reg = linregress(np.log(x), np.log(y))
    nearest_slope = find_nearest(np.arange(-7, 7, 0.25), reg.slope)
    if loc == "up":
        sign = 1
        valign = "bottom"
        halign = "left"
    else:
        sign = -1
        valign = "top"
        halign = "right"
    dx = x[-1] - x[-2]
    dy = y[-1] - y[-2]
    x1 = x[-2] + 0.25*dx
    x2 = x1 + 0.5*dx
    y1 = y[-2] + 0.25*dy + 0.2*sign*abs(dy)
    y2 = np.exp( nearest_slope*np.log(x2/x1) + np.log(y1) )
    ax.plot([x1,x2], [y1,y2], "k-")
    ax.text(
        np.average(0.5*(x1+x2) + 0.05*sign*abs(dx)),
        np.average(0.5*(y1+y2) + 0.05*sign*abs(dy)),
        r"${}$".format(abs(nearest_slope)),
        ha=halign,
        va=valign
    )



print("Doing case analysis in {}".format(os.getcwd()))

steps_to_do = [1,2,3]

# directory where outsourced scripts lie
script_dir = "/home/vachan/Documents/Work/plens/study1/riemann_1d/scripts/"

# constant data (generally not required to change this)
N_values = [1,2,3,5]
dofs = [200, 400, 800]
actual_dof_values = pd.DataFrame(
    [
        [100*2, 200*2, 400*2], # N=1
        [66*3, 134*3, 266*3], # N=2
        [50*4, 100*4, 200*4], # N=3
        [34*6, 66*6, 134*6] # N=5
    ],
    index=N_values,
    columns=dofs
)

# varying data (requires user intervention)
flux = "chandrashekhar"
flux_display = "Chandrashekhar" # how the flux scheme should be printed/written on plots
result_dir = "result_logarithm"
case_name = "Test 5"
individual_analysis_file = "full_analysis.log"



# 1. Group plots showing results with different N for dixed dof (outsourced)
if 1 in steps_to_do:
    for dof in dofs:
        subprocess.run([
            "python3",
            script_dir + "plot_group_dof.py",
            ".",
            str(dof),
            flux,
            result_dir,
            "comparison_data.csv",
            # "{}, {} dof, {}".format(case_name, dof, flux_display),
            "",
            "--save",
            "../plots",
            "dof{}_12_12_{}".format(dof, flux),
            "--size",
            "7",
            "6"
        ])



# 2. Error vs dof for different values of N
# Create a data frames for this purpose
l2_errors = pd.DataFrame(
    np.zeros((len(N_values), len(dofs))),
    index=N_values,
    columns=dofs
)
wtpt = l2_errors.copy() # wall time per time step
ctpt = l2_errors.copy() # cpu time per time step
for N in N_values:
    for dof in dofs:
        df = pd.read_csv(
            "dof{}_12_12_N{}_{}/{}/{}".format(
                dof, N, flux, result_dir, individual_analysis_file
            ),
            index_col=0,
            header=None
        )
        # clean the data: keep only latest logged data
        df = df[~df.index.duplicated(keep="last")]
        l2_errors.loc[N, dof] = df.loc["l2 error", 1]
        wtpt.loc[N, dof] = df.loc["wall time per time step", 1]
        ctpt.loc[N, dof] = df.loc["cpu time per time step", 1]

N_linestyles = pd.Series(["-", "--", "-.", ":"], index=N_values)
N_markers = pd.Series(["o", "s", "P", "^"], index=N_values)
N_markercolors = pd.Series(["blue", "green", "red", "magenta"], index=N_values)
dof_linestyles = pd.Series(["-", "--", "-."], index=dofs)

if 2 in steps_to_do:
    fig, ax = plt.subplots(1,1)
    for N in N_values:
        ax.plot(
            actual_dof_values.loc[N, :],
            l2_errors.loc[N, :],
            ls=N_linestyles.loc[N],
            marker=N_markers[N],
            c=N_markercolors[N],
            # markerfacecolor=N_markercolors[N],
            # markeredgecolor=N_markercolors[N],
            label=r"$N={}$".format(N)
        )
    plot_convergence_rate(ax, actual_dof_values.loc[5, :], l2_errors.loc[5, :], "below")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    # ax.xaxis.set_minor_locator(MultipleLocator(50))
    # ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_locator(MultipleLocator(2e-2))
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_locator(MultipleLocator(5e-3))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel("Degrees of freedom")
    ax.set_ylabel(r"$L^2$ error")
    # ax.set_title("{}, {}".format(case_name, flux_display))
    ax.legend(loc="best", handlelength=3)
    ax.grid(which="major")
    fig.set_size_inches(4, 4)
    fig.tight_layout(rect=[0,0,1,1])
    plt.show()
    mysavefig(fig, "../plots", "error_vs_dof_{}".format(flux))



# 3. Plot error vs cpu time for different N and across dofs
if 3 in steps_to_do:
    fig, ax = plt.subplots(1,1)
    for N in N_values:
        ax.plot(
            ctpt.loc[N, :],
            l2_errors.loc[N, :],
            ls=N_linestyles.loc[N],
            marker=N_markers[N],
            c=N_markercolors[N],
            # markerfacecolor=N_markercolors[N],
            # markeredgecolor=N_markercolors[N],
            label=r"$N={}$".format(N)
        )
    ax.set_xlabel("CPU time per time step [sec]")
    ax.set_ylabel(r"$L^2$ error")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_locator(MultipleLocator(2e-2))
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_locator(MultipleLocator(5e-3))
    ax.yaxis.set_minor_formatter(NullFormatter())
    # ax.set_title("{}, {}".format(case_name, flux_display))
    ax.legend(loc="best", handlelength=3)
    ax.grid(which="major")
    fig.set_size_inches(4, 4)
    fig.tight_layout(rect=[0,0,1,1])
    plt.show()
    mysavefig(fig, "../plots", "error_vs_cputime_{}".format(flux))



# 4. Visual comparison for all dofs and all values of N
# Don't look so good though
if 4 in steps_to_do:
    subprocess.run([
        "python3",
        script_dir + "plot_group_all.py",
        ".",
        flux,
        result_dir,
        "comparison_data.csv",
        "{}, {}".format(case_name, flux_display),
        "--save",
        "../plots",
        "all_{}".format(flux),
        "--size",
        "7",
        "6"
    ])