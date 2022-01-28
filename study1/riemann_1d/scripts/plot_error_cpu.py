# Plots errors in multiple cases with varying N gives a scatter plot of cpu time vs error

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter, NullFormatter, MultipleLocator, LogLocator
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 12
plt.rcParams["mathtext.fontset"] = "dejavuserif"

dof_values = [200, 400]

case_names = [
    "Test 1, HLLC",
    "Test 2, Rusanov",
    "Test 3, HLLC",
    "Test 4, HLLC",
    "Test 5, Rusanov",
]

N_values = [1,2,3,5]

# errors calculated from plot_individual.py
# ||simulation-exact|| / || exact ||
# access: l1_errors[<dof>][<case>][<N>]
l1_errors = np.array([
    [ # 200 dofs
        [1.693e-2, 1.078e-2, 0.7878e-2, 0.7375e-2], # case 1
        [0.2311, 0.1655, 0.1235, 0.09526],
        [0.1565, 0.1253, 0.1232, 0.1003],
        [4.391e-2, 4.053e-2, 3.582e-2, 3.69e-2],
        [0.1399, 0.09716, 0.1027, 0.0956] # case 5
    ],

    [ # 400 dofs
        [1.0e-2, 0.598e-2, 0.424e-2, 0.4755e-2],
        [0.1439, 0.09731, 0.07329, 0.0576],
        [0.09814, 0.06594, 0.06426, 0.05333],
        [2.811e-2, 2.553e-2, 2.225e-2, 2.518e-2],
        [8.44e-2, 5.228e-2, 5e-2, 4.725e-2]
    ]
])
l2_errors = np.array([
    [
        [3.55e-2, 2.867e-2, 2.446e-2, 2.225e-2],
        [0.2858, 0.2441, 0.1819, 0.139],
        [0.3609, 0.3322, 0.3165, 0.2766],
        [11.84e-2, 12.64e-2, 11.01e-2, 10.34e-2],
        [0.3314, 0.2702, 0.2735, 0.2633]
    ],

    [
        [2.825e-2, 2.114e-2, 1.708e-2, 1.581e-2],
        [0.2026, 0.1733, 0.1305, 0.1013],
        [0.2829, 0.2386, 0.2255, 0.2],
        [9.196e-2, 9.255e-2, 8.845e-2, 7.404e-2],
        [25.91e-2, 19.94e-2, 19.78e-2, 17.16e-2]
    ]
])
linf_errors = np.array([
    [
        [0.1545, 0.1337, 0.1384, 0.1413],
        [0.4605, 0.4963, 0.4353, 0.3991],
        [0.61, 0.6171, 0.6441, 0.6057],
        [60.42e-2, 60.02e-2, 57.35e-2, 53.8e-2],
        [0.6261, 0.5716, 0.5974, 0.6606]
    ],

    [
        [15.12e-2, 13.91e-2, 14.28e-2, 13.56e-2],
        [0.4243, 0.4565, 0.3978, 0.3794],
        [0.6691, 0.6981, 0.6812, 0.624],
        [59.86e-2, 57.19e-2, 70.57e-2, 54.07e-2],
        [60.7e-2, 58.39e-2, 74.1e-2, 61.47e-2]
    ]
])
all_errors = [l1_errors, l2_errors, linf_errors]
all_errors_types = ["L1", "L2", "Linf"]

# cpu time per time step
cpu_times = np.array([
    [
        [9.847, 8.917, 9.346, 14.23],
        [10.48, 9.391, 9.793, 13.37],
        [11.29, 9.736, 9.811, 14.8],
        [10.73, 9.648, 10.73, 15.07],
        [17.03, 9.287, 8.914, 13.43]
    ],

    [
        [21.71, 17.38, 17.82, 27.36],
        [19.86, 17.79, 19.78, 27.75],
        [21.14, 21.67, 17.92, 25.51],
        [19.72, 18.01, 18.26, 25.47],
        [19.44, 23.17, 19.35, 26.31]
    ]
])



# Error vs N
"""
case_markers = ["o", "s", "^", "d", "x"]
case_colors = ["b", "r", "g", "m", "k"]
dof_linestyles = ["-", "--"]
for error_type, errors in zip(all_errors_types, all_errors):
    fig, ax = plt.subplots(1,1)
    for dof_id in range(len(dof_values)):
        for case_id in range(len(case_names)):
            # plot error scaled with N=1 value
            ax.plot(
                N_values,
                np.array(errors[dof_id,case_id])/errors[dof_id,case_id,0],
                marker=case_markers[case_id],
                c=case_colors[case_id],
                ls=dof_linestyles[dof_id]
            )
    # loop for legend
    legend_elements = []
    for case_id in range(len(case_names)):
        legend_elements.append(
            Line2D(
                [0], [0],
                marker=case_markers[case_id],
                markerfacecolor=case_colors[case_id],
                markeredgecolor=case_colors[case_id],
                ls="",
                label=case_names[case_id]
            )
        )
    for dof_id in range(len(dof_values)):
        legend_elements.append(
            Line2D(
                [0], [0],
                c="gray",
                ls=dof_linestyles[dof_id],
                label="{} dofs".format(dof_values[dof_id])
            )
        )
    ax.grid()
    ax.legend(handles=legend_elements, loc="best")
    ax.set_xlabel(r"$N$")
    ax.set_ylabel("Scaled {} error".format(error_type))
    fig.tight_layout(rect=[0,0,1,1])
    def on_resize(event):
        fig.tight_layout(rect=[0,0,1,1])
        fig.canvas.draw()
    cid = fig.canvas.mpl_connect('resize_event', on_resize)
    plt.show()
    del fig, ax
"""



# CPU time vs error
# https://stackoverflow.com/questions/12761806/matplotlib-2-different-legends-on-same-graph
N_markers = ["o", "s", "P", "^"]
N_markercolors = ["blue", "green", "red", "magenta"]
case_connector_linestyles = [
    (0, ()), # solid line
    (0, (5,5)), # regular dash
    (0, (5,1)), # dense dash
    (0, (5,2,1,2)), # dash dot
    (0, (1,1)) # dot
]
dof_linecolors = ["orange", "brown"]

# first draw the connectors
fig, ax = plt.subplots(1,1)
for dof_id in range(len(dof_values)):
    for case_id in range(len(case_names)):
        ax.plot(
            np.array(l2_errors[dof_id,case_id,:])/l2_errors[0,case_id,0], # scale by N=1 error of 200 dof error
            np.array(cpu_times[dof_id,case_id,:]),#/cpu_times[dof_id][case_id][0], # scale by N=1 cpu time
            c=dof_linecolors[dof_id],
            ls=case_connector_linestyles[case_id],
        )
# now add the markers
for N_id in range(len(N_values)):
    for dof_id in range(len(dof_values)):
        ax.plot(
            np.array(l2_errors[dof_id,:,N_id])/np.array(l2_errors[0,:,0]),
            np.array(cpu_times[dof_id,:,N_id]),
            ls="",
            marker=N_markers[N_id],
            markerfacecolor=N_markercolors[N_id],
            markeredgecolor=N_markercolors[N_id],
        )
# now make the legend
legend_elements = []
for N_id in range(len(N_values)):
    legend_elements.append(
        Line2D(
            [0], [0],
            ls="",
            marker=N_markers[N_id],
            markerfacecolor=N_markercolors[N_id],
            markeredgecolor=N_markercolors[N_id],
            label=r"$N={}$".format(N_values[N_id])
        )
    )
for case_id in range(len(case_names)):
    legend_elements.append(
        Line2D(
            [0], [0],
            ls=case_connector_linestyles[case_id],
            c="gray",
            label=case_names[case_id]
        )
    )
for dof_id in range(len(dof_values)):
    legend_elements.append(
        Line2D(
            [0], [0],
            c=dof_linecolors[dof_id],
            label="{} dofs".format(dof_values[dof_id])
        )
    )

ax.legend(handles=legend_elements, loc="upper right", bbox_to_anchor=(1,1))
ax.set_xlabel("Scaled normalised L2 error")
ax.set_ylabel("CPU time per time step [sec]")
ax.set_xscale("log")
ax.set_yscale("log")
ax.xaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.xaxis.set_minor_formatter(NullFormatter())
ax.yaxis.set_minor_formatter(NullFormatter())
ax.grid(which="major", lw=0.5)
ax.grid(which="minor", lw=0.25)
fig.tight_layout()
def on_resize(event):
    fig.tight_layout()
    fig.canvas.draw()
cid = fig.canvas.mpl_connect('resize_event', on_resize)
plt.show()
