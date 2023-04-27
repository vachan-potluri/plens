# A script to do a complete analysis of an individual test case where multiple runs with varying
# dof and N have been done. It is assumed that for a given flux scheme, a case has been run with
# 200, 400 and 800 dofs in x-dir (nominally) and 12 dofs in y and z directions. It is also assumed
# that the N values are 1, 2, 3 and 5. The script generates the following plots
#
# 1. Visual comparison of results for given dof and varying N (outsourced)
# 2. Error vs dof for different values of N
# 3. Error vs cpu time per time step for different values of N and across dofs
# 4. Visual comparison of results for all dofs and all values of N (outsourced)
# 5. Error vs N for different values of dof (was suggested by KB, see WJ-04-Mar-2022)
# 51. Error*DoF vs N (for PMRF review May 2023)
# 6. CPU time/DoF/time step vs N (for SIS2022, see WJ-16-Jun2022)
# 7. CPU time/DoF vs N (for APS-3 presentation, see WJ-20-Sep-2022)
# 8. Error vs cpu time for different values of N and across dofs (for APS-3 presentation, see WJ-20-Sep-2022)
# 9. Total variation vs N for different values of dof (see WJ-21-Mar-2023)
# 91. A second formula for total variation (see WJ-23-Mar-2023)
# 10. Variations of cumulative total variations (see WJ-22-Mar-2023)
# 11. Plot error vs wall time for different values of N across dofs (for manuscript2, see plot_error_cpu_vs.py)
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
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter, NullFormatter, MultipleLocator
from scipy.stats import linregress
plt.rcParams.update({
    "axes.formatter.limits": [-2,2],
    "font.family": "Times",
    "font.size": 12,
    "text.latex.preamble": r"\usepackage{siunitx}",
    "text.usetex": True,
})



def calc_total_variation(array):
    # calculates total variation of an array
    tv = 0.0
    for i in range(len(array)-1):
        tv += abs(array[i+1] - array[i])
    return tv



def get_cumulative_total_variation(array):
    # calculates cumulative total variation: total variation of a segment from the start to current
    # location
    tv_cum = np.zeros_like(array)
    for i in range(len(array)-1):
        tv_cum[i+1] = tv_cum[i] + abs(array[i+1] - array[i])
    return tv_cum



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



def format_error_axis(axis, major_loc, minor_loc):
    # formats the error axis which is common for steps 2, 3 and 5
    # ax.yaxis can be passed as the argument to this function with the major and minor tick
    # separations provided as arguments to this function
    axis.set_major_locator(MultipleLocator(major_loc))
    axis.set_major_formatter(ScalarFormatter())
    axis.set_minor_locator(MultipleLocator(minor_loc))
    axis.set_minor_formatter(NullFormatter())



def format_dof_axis(axis, major_loc=100, minor_loc=None):
    # similar to `format_error_axis`, but for the dof axis
    # minor ticks generally not required here, but are added if the parameter passed is not None
    ax.xaxis.set_major_locator(MultipleLocator(major_loc))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    if minor_loc is not None:
        ax.xaxis.set_minor_locator(MultipleLocator(minor_loc))
        ax.xaxis.set_minor_formatter(NullFormatter())


def format_ctpt_axis(axis, major_loc=1.0, minor_loc=0.25):
    # similar to `format_error_axis`, but for the cpu time/time step axis
    axis.set_major_locator(MultipleLocator(major_loc))
    axis.set_major_formatter(ScalarFormatter())
    axis.set_minor_locator(MultipleLocator(minor_loc))
    axis.set_minor_formatter(NullFormatter())



print("Doing case analysis in {}".format(os.getcwd()))

steps_to_do = [51]

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

N_linestyles = pd.Series(["-", "--", "-.", ":"], index=N_values)
N_markers = pd.Series(["o", "s", "P", "^"], index=N_values)
N_markercolors = pd.Series(["blue", "green", "red", "magenta"], index=N_values)
dof_linestyles = pd.Series(["-", "--", "-."], index=dofs)
dof_linecolors = pd.Series(["r", "g", "b"], index=dofs)
dof_markers = pd.Series(["o", "s", "^"], index=dofs)

# varying data (requires user intervention)
flux = "chandrashekhar"
flux_display = "Chandrashekar" # how the flux scheme should be printed/written on plots
result_dir = "result_16Jan2023_persson"
# result_dir = "result_28Mar2022"
# result_dir = "result_04Apr2023_timing"
case_name = "Test 1"
individual_analysis_file = "full_analysis.log"
# major and minor locators for error axis, changes on test-by-test basis
error_ax_major_loc = 5e-3 # test1-1
# error_ax_major_loc = 1e-1 # test1-2
# error_ax_major_loc = 5e-2 # test1-3
# error_ax_major_loc = 1e-2 # test1-4
# error_ax_major_loc = 3e-2 # test1-5
error_ax_minor_loc = error_ax_major_loc/4



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
            "5",
            "3"
        ])



# 2. Error vs dof for different values of N
# Create a data frames for this purpose
l1_errors = pd.DataFrame(
    np.zeros((len(N_values), len(dofs))),
    index=N_values,
    columns=dofs
)
wtpt = l1_errors.copy() # wall time per time step
ctpt = l1_errors.copy() # cpu time per time step
cpu_time = l1_errors.copy() # cpu time
wall_times = l1_errors.copy()
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
        l1_errors.loc[N, dof] = df.loc["l1 error", 1]
        wtpt.loc[N, dof] = df.loc["wall time per time step", 1]
        ctpt.loc[N, dof] = df.loc["cpu time per time step", 1]
        cpu_time.loc[N, dof] = df.loc["cpu time", 1]
        wall_times.loc[N, dof] = wtpt.loc[N, dof]*df.loc["time steps", 1]

if 2 in steps_to_do:
    fig, ax = plt.subplots(1,1)
    for N in N_values:
        ax.plot(
            actual_dof_values.loc[N, :],
            l1_errors.loc[N, :],
            ls=N_linestyles.loc[N],
            marker=N_markers[N],
            c=N_markercolors[N],
            # markerfacecolor=N_markercolors[N],
            # markeredgecolor=N_markercolors[N],
            label=r"$N={}$".format(N)
        )
    plot_convergence_rate(ax, actual_dof_values.loc[5, :], l1_errors.loc[5, :], "below")
    ax.set_xscale("log")
    ax.set_yscale("log")
    format_dof_axis(ax.xaxis)
    format_error_axis(ax.yaxis, error_ax_major_loc, error_ax_minor_loc)
    ax.set_xlabel("Degrees of freedom")
    ax.set_ylabel(r"$L^1$ error")
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
            l1_errors.loc[N, :],
            ls=N_linestyles.loc[N],
            marker=N_markers[N],
            c=N_markercolors[N],
            # markerfacecolor=N_markercolors[N],
            # markeredgecolor=N_markercolors[N],
            label=r"$N={}$".format(N)
        )
    ax.set_xlabel("CPU time/time step [sec]")
    ax.set_ylabel(r"$L^1$ error")
    ax.set_xscale("log")
    ax.set_yscale("log")
    format_ctpt_axis(ax.xaxis)
    format_error_axis(ax.yaxis, error_ax_major_loc, error_ax_minor_loc)
    # ax.set_title("{}, {}".format(case_name, flux_display))
    ax.legend(loc="best", handlelength=3)
    ax.grid(which="major")
    fig.set_size_inches(5, 3.5)
    fig.tight_layout(rect=[0,0,1,1])
    plt.show()
    mysavefig(fig, "../plots", "error_vs_cputime_per_timestep_{}".format(flux))



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



# 5. Error vs N for different dof values
if 5 in steps_to_do:
    fig, ax = plt.subplots(1,1)
    for dof in dofs:
        ax.plot(
            N_values,
            l1_errors.loc[:, dof],
            ls=dof_linestyles.loc[dof],
            marker=dof_markers.loc[dof],
            c=dof_linecolors.loc[dof],
            # markerfacecolor=N_markercolors[N],
            # markeredgecolor=N_markercolors[N],
            label=r"{} DoFs".format(dof)
        )
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"Normalised $L^1$ error")
    ax.set_yscale("log")
    format_error_axis(ax.yaxis, error_ax_major_loc, error_ax_minor_loc)
    # ax.set_title("{}, {}".format(case_name, flux_display))
    ax.legend(loc="best", handlelength=3)
    ax.grid(which="major")
    ax.xaxis.set_major_locator(MultipleLocator(1)) # don't need non-integer values for N
    fig.set_size_inches(5, 2.5)
    fig.tight_layout(rect=[0,0,1,1], pad=0.25)
    plt.show()
    mysavefig(fig, "../plots", "error_vs_N_{}".format(flux))



# 51. Error*DoFs vs N
if 51 in steps_to_do:
    fig, ax = plt.subplots(1,1)
    for dof in dofs:
        ax.plot(
            N_values,
            l1_errors.loc[:, dof]*dof,
            ls=dof_linestyles.loc[dof],
            marker=dof_markers.loc[dof],
            c=dof_linecolors.loc[dof],
            # markerfacecolor=N_markercolors[N],
            # markeredgecolor=N_markercolors[N],
            label=r"{} DoFs".format(dof)
        )
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"Error*DoFs")
    ax.set_yscale("log")
    ax.legend(loc="best", handlelength=3)
    ax.grid(which="major")
    ax.xaxis.set_major_locator(MultipleLocator(1)) # don't need non-integer values for N
    fig.set_size_inches(5, 2.5)
    fig.tight_layout(rect=[0,0,1,1], pad=0.25)
    plt.show()
    # mysavefig(fig, "../plots", "error_per_dof_vs_N_{}".format(flux))




# 6. CPU time/DoF/time step vs N
if 6 in steps_to_do:
    fig, ax = plt.subplots(1,1)
    # plot
    for N in N_values:
        for dof in dofs:
            ax.scatter(
                N,
                ctpt.loc[N, dof]/dof,
                marker=dof_markers.loc[dof],
                edgecolor=dof_linecolors.loc[dof],
                facecolor="none"
            )
    # for legend
    legend_elements = []
    for dof in dofs:
        legend_elements.append(
            Line2D(
                [0], [0],
                ls="",
                marker=dof_markers.loc[dof],
                markeredgecolor=dof_linecolors.loc[dof],
                markerfacecolor="none",
                label="{} DoFs".format(dof)
            )
        )
    ax.set_xlabel(r"$N$")
    ax.set_ylabel("CPU time/time step/DoF [sec]")
    ax.grid()
    ax.legend(handles=legend_elements)
    ax.xaxis.set_major_locator(MultipleLocator(1)) # don't need non-integer values for N
    fig.set_size_inches(4, 3)
    fig.tight_layout(rect=[0,0,1,1], pad=0.25)
    plt.show()
    mysavefig(fig, "../plots", "cputime_per_timestep_vs_N_{}".format(flux))



# 7. CPU time/DoF vs N (for APS-3 presentation, see WJ-20-Sep-2022)
if 7 in steps_to_do:
    fig, ax = plt.subplots(1,1)
    # plot
    for N in N_values:
        for dof in dofs:
            ax.scatter(
                N,
                cpu_time.loc[N, dof]/cpu_time.loc[1, dof],
                marker=dof_markers.loc[dof],
                edgecolor=dof_linecolors.loc[dof],
                facecolor="none"
            )
    # for legend
    legend_elements = []
    for dof in dofs:
        legend_elements.append(
            Line2D(
                [0], [0],
                ls="",
                marker=dof_markers.loc[dof],
                markeredgecolor=dof_linecolors.loc[dof],
                markerfacecolor="none",
                label="{} DoFs".format(dof)
            )
        )
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"CPU time($N, \bullet$)/CPU time($1, \bullet$)")
    ax.grid()
    ax.legend(handles=legend_elements)
    ax.xaxis.set_major_locator(MultipleLocator(1)) # don't need non-integer values for N
    fig.set_size_inches(4, 3)
    fig.tight_layout(rect=[0,0,1,1], pad=0.25)
    plt.show()
    mysavefig(fig, "../plots", "cputime_vs_N_{}".format(flux))



# 8. Error vs cpu time for different values of N and across dofs (for APS-3 presentation, see WJ-20-Sep-2022)
if 8 in steps_to_do:
    fig, ax = plt.subplots(1,1)
    # legends
    legend_elements = []
    # draw connecting lines for the same dof count
    for dof in dofs:
        line, = ax.plot(
            cpu_time.loc[:, dof]/cpu_time.loc[1, dof],
            l1_errors.loc[:, dof]/l1_errors.loc[1, dof],
            ls=dof_linestyles.loc[dof],
            c="gray",
            label="{} DoFs".format(dof)
        )
        legend_elements.append(line)
    # y=1/x curve
    x = np.linspace(1, 2)
    line, = ax.plot(x, 1/x, c="brown", ls=":", label=r"$y=1/x$")
    legend_elements.append(line)
    # add data points
    for N in N_values:
        line, = ax.plot(
            cpu_time.loc[N, :]/cpu_time.loc[1, :],
            l1_errors.loc[N, :]/l1_errors.loc[1, :],
            ls="",
            marker=N_markers[N],
            c=N_markercolors[N],
            label=r"$N={}$".format(N)
        )
        legend_elements.append(line)
    # for N in N_values:
    #     legend_elements.append(
    #         Line2D(
    #             [0], [0],
    #             ls="",
    #             marker=N_markers[N],
    #             c=N_markercolors[N],
    #             label=r"$N={}$".format(N)
    #         )
    #     )
    # for dof in dofs:
    #     legend_elements.append(
    #         Line2D(
    #             [0], [0],
    #             ls=dof_linestyles.loc[dof],
    #             c="gray",
    #             label="{} DoFs".format(dof)
    #         )
    #     )
    # legend_elements.append(
    #     Line2D(
    #         [0], [0], c="brown", ls=":", label=r"$y=1/x$"
    #     )
    # )
    ax.set_xlabel(r"CPU time($N, \bullet$)/CPU time($1, \bullet$)")
    ax.set_ylabel(r"$L^1$ error($N, \bullet$)/$L^1$ error($1, \bullet$)")
    ax.legend(handles=legend_elements, ncol=2)
    ax.grid(which="major")
    fig.set_size_inches(4, 3)
    fig.tight_layout(rect=[0,0,1,1], pad=0.25)
    plt.show()
    mysavefig(fig, "../plots", "error_vs_cputime_{}".format(flux))



# 9. Calculate and plot the total variation of the simulation error
if 9 in steps_to_do:
    tv = pd.DataFrame(index=N_values, columns=dofs)
    for dof in dofs:
        for N in N_values:
            data = np.genfromtxt(
                f"dof{dof}_12_12_N{N}_{flux}/{result_dir}/comparison_data.csv",
                delimiter=","
            )
            tv_sim = calc_total_variation(data[:,2])
            tv_ex = calc_total_variation(data[:,1])
            tv.loc[N, dof] = (tv_sim-tv_ex)/tv_ex
    fig, ax = plt.subplots(1,1)
    for dof in dofs:
        ax.plot(
            N_values,
            tv.loc[:, dof],
            ls=dof_linestyles.loc[dof],
            marker=dof_markers.loc[dof],
            c=dof_linecolors.loc[dof],
            # markerfacecolor=N_markercolors[N],
            # markeredgecolor=N_markercolors[N],
            label=r"{} DoFs".format(dof)
        )
    ax.set_xlabel(r"$N$")
    ax.set_ylabel("Normalised total variation")
    ax.legend(loc="best", handlelength=3)
    ax.grid()
    ax.xaxis.set_major_locator(MultipleLocator(1)) # don't need non-integer values for N
    fig.set_size_inches(5, 2.5)
    fig.tight_layout(rect=[0,0,1,1], pad=0.25)
    plt.show()
    mysavefig(fig, "../plots", "tv_vs_N_{}".format(flux))



# 91. Calculate and plot the total variation of the simulation error
if 91 in steps_to_do:
    tv = pd.DataFrame(index=N_values, columns=dofs)
    tv2 = tv.copy()
    for dof in dofs:
        for N in N_values:
            data = np.genfromtxt(
                f"dof{dof}_12_12_N{N}_{flux}/{result_dir}/comparison_data.csv",
                delimiter=","
            )
            tv_ex = calc_total_variation(data[:,1])
            tv_sim = calc_total_variation(data[:,2])
            tv.loc[N, dof] = (tv_sim-tv_ex)/tv_ex
            tv2.loc[N, dof] = calc_total_variation(data[:,2]-data[:,1])
    fig, axes = plt.subplots(1,2)
    for dof in dofs:
        axes[0].plot(
            N_values,
            tv.loc[:, dof],
            ls=dof_linestyles.loc[dof],
            marker=dof_markers.loc[dof],
            c=dof_linecolors.loc[dof],
            # markerfacecolor=N_markercolors[N],
            # markeredgecolor=N_markercolors[N],
            label=r"{} DoFs".format(dof)
        )
        axes[0].set_ylabel("Normalised TV difference")
        axes[1].plot(
            N_values,
            tv2.loc[:, dof],
            ls=dof_linestyles.loc[dof],
            marker=dof_markers.loc[dof],
            c=dof_linecolors.loc[dof],
            # markerfacecolor=N_markercolors[N],
            # markeredgecolor=N_markercolors[N],
            label=r"{} DoFs".format(dof)
        )
        axes[1].set_ylabel("TV of the error")
    for ax in axes:
        ax.set_xlabel(r"$N$")
        ax.legend(loc="best", handlelength=3)
        ax.grid()
        ax.xaxis.set_major_locator(MultipleLocator(1)) # don't need non-integer values for N
    fig.set_size_inches(6, 2.5)
    fig.tight_layout(rect=[0,0,1,1], pad=0.25)
    plt.show()
    mysavefig(fig, "../plots", "tv_tv2_vs_N_{}".format(flux))



# 10: cumulative total variation
if 10 in steps_to_do:
    fig, ax = plt.subplots(1,1,figsize=(8,6))
    for dof in [200,400,800]:
        for N in N_values:
            data = np.genfromtxt(
                f"dof{dof}_12_12_N{N}_{flux}/{result_dir}/comparison_data.csv",
                delimiter=","
            )
            tv_cum = get_cumulative_total_variation(data[:,2])
            ax.plot(data[:,0], tv_cum, label=rf"{dof} DoF, $N$={N}")
    ax.set_xlabel("$x$")
    ax.set_ylabel(r"$\textrm{TV}(\phi(\xi))$ for $\xi \in [0,x]$")
    ax.grid()
    ax.legend()
    fig.tight_layout()
    plt.show()



# 11. Error vs wall time for manuscript2
if 11 in steps_to_do:
    fig, ax = plt.subplots(1,1,figsize=(6,4))
    for N in N_values:
        ax.plot(
            wall_times.loc[N,:],
            l1_errors.loc[N,:],
            ls=N_linestyles.loc[N],
            c=N_markercolors.loc[N],
            marker=N_markers.loc[N],
            label=rf"$N={N}$"
        )
    ax.grid(which="major", lw=1)
    ax.grid(which="minor", lw=0.5)
    ax.legend(handlelength=3)
    ax.set_xlabel(r"Wall-clock time [\unit{s}]")
    ax.set_ylabel(r"Normalised $L^1$ error")
    ax.set_xscale("log")
    ax.set_yscale("log")
    fig.tight_layout(pad=0.25)
    plt.show()
