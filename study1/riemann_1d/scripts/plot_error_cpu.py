# Plots errors in multiple cases with varying N gives a scatter plot of cpu time vs error

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 10
plt.rcParams["mathtext.fontset"] = "dejavuserif"

case_names = np.array([
    "Test 1, 200 dof, HLLC",
    "Test 2, 200 dof, Rusanov",
    "Test 3, 200 dof, HLLC",
    "Test 4*, 200 dof, HLLC",
    "Test 5, 200 dof, Rusanov",
    
    "Test 1*, 400 dof, HLLC",
    "Test 2, 400 dof, Rusanov",
    "Test 3, 400 dof, HLLC",
    "Test 4*, 400 dof, HLLC"
])

case_legend_names = np.array([
    "Test 1, 200 dof, HLLC",
    "Test 2, 200 dof, Rusanov",
    "Test 3, 200 dof, HLLC",
    "Test 4*, 200 dof, HLLC",
    "Test 5, 200 dof, Rusanov",

    "Test 1*, 400 dof, HLLC",
    "Test 2, 400 dof, Rusanov",
    "Test 3, 400 dof, HLLC",
    "Test 4*, 400 dof, HLLC"
])

# errors calculated from plot_individual.py
# ||simulation-exact|| / || exact ||
N_values = np.array([1,2,3,5])
# access: l1_errors[<case>][<N>]
l1_errors = np.array([
    [1.693e-2, 1.078e-2, 0.7878e-2, 0.7375e-2], # first case in case_names: errors corresponding to N_values
    [0.2311, 0.1655, 0.1235, 0.09526],
    [0.1565, 0.1253, 0.1232, 0.1003],
    [4.391e-2, 4.053e-2, 3.582e-2, 3.69e-2],
    [0.1399, 0.09716, 0.1027, 0.0956],

    [1.0e-2, 0.598e-2, 0.424e-2, 0.4755e-2],
    [0.1439, 0.09731, 0.07329, 0.0576],
    [0.09814, 0.06594, 0.06426, 0.05333],
    [2.811e-2, 2.553e-2, 2.225e-2, 2.518e-2] # last case in case_names
])
l2_errors = np.array([
    [3.55e-2, 2.867e-2, 2.446e-2, 2.225e-2],
    [0.2858, 0.2441, 0.1819, 0.139],
    [0.3609, 0.3322, 0.3165, 0.2766],
    [11.84e-2, 12.64e-2, 11.01e-2, 10.34e-2],
    [0.3314, 0.2702, 0.2735, 0.2633],

    [2.825e-2, 2.114e-2, 1.708e-2, 1.581e-2],
    [0.2026, 0.1733, 0.1305, 0.1013],
    [0.2829, 0.2386, 0.2255, 0.2],
    [9.196e-2, 9.255e-2, 8.845e-2, 7.404e-2]
])
linf_errors = np.array([
    [0.1545, 0.1337, 0.1384, 0.1413],
    [0.4605, 0.4963, 0.4353, 0.3991],
    [0.61, 0.6171, 0.6441, 0.6057],
    [60.42e-2, 60.02e-2, 57.35e-2, 53.8e-2],
    [0.6261, 0.5716, 0.5974, 0.6606],

    [15.12e-2, 13.91e-2, 14.28e-2, 13.56e-2],
    [0.4243, 0.4565, 0.3978, 0.3794],
    [0.6691, 0.6981, 0.6812, 0.624],
    [59.86e-2, 57.19e-2, 70.57e-2, 54.07e-2]
])
all_errors = np.array([l1_errors, l2_errors, linf_errors])
all_errors_types = np.array(["L1", "L2", "Linf"])

# cpu time per time step
cpu_times = np.array([
    [9.847, 8.917, 9.346, 14.23],
    [10.48, 9.391, 9.793, 13.37],
    [11.29, 9.736, 9.811, 14.8],
    [10.73, 9.648, 10.73, 15.07],
    [17.03, 9.287, 8.914, 13.43],

    [21.71, 17.38, 17.82, 27.36],
    [19.86, 17.79, 19.78, 27.75],
    [21.14, 21.67, 17.92, 25.51],
    [19.72, 18.01, 18.26, 25.47]
])



# Error vs N
"""
case_plot_formats = np.array([
    "bo-",
    "rs-",
    "g^-",
    "md-",
    "kx-"
])
for error_type, errors in zip(all_errors_types, all_errors):
    fig, ax = plt.subplots(1,1)
    for case_id in range(len(case_names)):
        # plot error scaled with N=1 value
        ax.plot(
            N_values,
            errors[case_id]/errors[case_id][0],
            case_plot_formats[case_id],
            label=case_legend_names[case_id]
        )
    ax.grid()
    ax.legend(loc="best")
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
N_markers = np.array(["bo", "gs", "rx", "m^"])
fig, ax = plt.subplots(1,1)
# connect data points related to a single case
connectors = []
connector_labels = case_legend_names
case_connector_colors = [
    "black", "black", "black", "black", "black",
    "brown", "brown", "brown", "brown", "brown"
]
case_connector_linestyles = [
    (0, ()), (0, (5,5)), (0, (5,1)), (0, (5,2,1,2)), (0, (1,1)),
    (0, ()), (0, (5,5)), (0, (5,1)), (0, (5,2,1,2)), (0, (1,1))
]
case_connector_linewidths = [
    1,1,1,1,1,
    1,1,1,1,1
]
for case_id in range(len(case_names)):
    l, = ax.plot(
        l2_errors[case_id,:]/l2_errors[case_id,0], # scale by N=1 error
        cpu_times[case_id,:],#/cpu_times[case_id,0], # scale by N=1 cpu time
        c=case_connector_colors[case_id],
        ls=case_connector_linestyles[case_id],
        lw=case_connector_linewidths[case_id],
        alpha = 0.75
    )
    connectors.append(l)
legend1 = ax.legend(connectors, connector_labels, loc="lower left", bbox_to_anchor=(1,0.5), handlelength=3)
ax.add_artist(legend1)
# plot individual data points
scatters = []
scatter_labels = ["N=1", "N=2", "N=3", "N=5"]
for N_id in range(len(N_values)):
    l, = ax.plot(
        l2_errors[:,N_id]/l2_errors[:,0], # scale by N=1 error
        cpu_times[:,N_id],#/cpu_times[:,0], # scale by N=1 cpu time
        N_markers[N_id]
    )
    scatters.append(l)
legend2 = ax.legend(scatters, scatter_labels, loc="upper left", bbox_to_anchor=(1,0.5))
# ax.add_artist(legend2)
ax.set_xlabel("Scaled L2 error")
ax.set_ylabel("CPU time per time step [sec]")
ax.grid()
fig.tight_layout()
def on_resize(event):
    fig.tight_layout()
    fig.canvas.draw()
cid = fig.canvas.mpl_connect('resize_event', on_resize)
fig.savefig('temp.pdf', bbox_extra_artists=(legend1,), bbox_inches='tight')
plt.show()
