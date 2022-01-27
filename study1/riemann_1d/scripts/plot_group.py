# Plots group plots using saved data from plot_individual.py

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 10
plt.rcParams["mathtext.fontset"] = "dejavuserif"

# data_files = [
#     "../test1-1/run/dof200_12_12_N1_hllc/result/comparison_data.csv",
#     "../test1-1/run/dof200_12_12_N2_hllc/result/comparison_data.csv",
#     "../test1-1/run/dof200_12_12_N3_hllc/result/comparison_data.csv",
#     "../test1-1/run/dof200_12_12_N5_hllc/result/comparison_data.csv"
# ]
# data_files = [
#     "../test1-2/run/dof200_12_12_N1_rusanov/result/comparison_data.csv",
#     "../test1-2/run/dof200_12_12_N2_rusanov/result/comparison_data.csv",
#     "../test1-2/run/dof200_12_12_N3_rusanov/result/comparison_data.csv",
#     "../test1-2/run/dof200_12_12_N5_rusanov/result/comparison_data.csv"
# ]
# data_files = [
#     "../test1-3/run/dof200_12_12_N1_hllc/result/comparison_data.csv",
#     "../test1-3/run/dof200_12_12_N2_hllc/result/comparison_data.csv",
#     "../test1-3/run/dof200_12_12_N3_hllc/result/comparison_data.csv",
#     "../test1-3/run/dof200_12_12_N5_hllc/result/comparison_data.csv"
# ]
# data_files = [
#     "../test1-4/run/dof200_12_12_N1_hllc/result/comparison_data.csv",
#     "../test1-4/run/dof200_12_12_N2_hllc/result/comparison_data.csv",
#     "../test1-4/run/dof200_12_12_N3_hllc/result/comparison_data.csv",
#     "../test1-4/run/dof200_12_12_N5_hllc/result/comparison_data.csv"
# ]
# data_files = [
#     "../test1-5/run/dof200_12_12_N1_rusanov/result/comparison_data.csv",
#     "../test1-5/run/dof200_12_12_N2_rusanov/result/comparison_data.csv",
#     "../test1-5/run/dof200_12_12_N3_rusanov/result/comparison_data.csv",
#     "../test1-5/run/dof200_12_12_N5_rusanov/result/comparison_data.csv"
# ]
# data_files = [
#     "../test1-1/run/dof400_12_12_N1_hllc/result/comparison_data.csv",
#     "../test1-1/run/dof400_12_12_N2_hllc/result/comparison_data.csv",
#     "../test1-1/run/dof400_12_12_N3_hllc/result/comparison_data.csv",
#     "../test1-1/run/dof400_12_12_N5_hllc/result1/comparison_data.csv"
# ]
# data_files = [
#     "../test1-2/run/dof400_12_12_N1_rusanov/result/comparison_data.csv",
#     "../test1-2/run/dof400_12_12_N2_rusanov/result/comparison_data.csv",
#     "../test1-2/run/dof400_12_12_N3_rusanov/result/comparison_data.csv",
#     "../test1-2/run/dof400_12_12_N5_rusanov/result/comparison_data.csv"
# ]
# data_files = [
#     "../test1-3/run/dof400_12_12_N1_hllc/result/comparison_data.csv",
#     "../test1-3/run/dof400_12_12_N2_hllc/result/comparison_data.csv",
#     "../test1-3/run/dof400_12_12_N3_hllc/result/comparison_data.csv",
#     "../test1-3/run/dof400_12_12_N5_hllc/result/comparison_data.csv"
# ]
data_files = [
    "../test1-4/run/dof400_12_12_N1_hllc/result/comparison_data.csv",
    "../test1-4/run/dof400_12_12_N2_hllc/result/comparison_data.csv",
    "../test1-4/run/dof400_12_12_N3_hllc/result/comparison_data.csv",
    "../test1-4/run/dof400_12_12_N5_hllc/result/comparison_data.csv"
]
titles = ["N=1", "N=2", "N=3", "N=5"]
figtitle = r"test1-4 $\rho$ vs $x$"

n_plots = len(data_files)
fig, axes = plt.subplots(2,2)
i = 0
for row in axes:
    for ax in row:
        data = np.genfromtxt(data_files[i], delimiter=",")
        ax.plot(data[:,0], data[:,1], "b-", label="Exact")
        ax.plot(data[:,0], data[:,2], "r-", label="Simulation")
        ax.set_title(titles[i])
        ax.grid()
        ax.legend(loc="best")
        i += 1
fig.suptitle(figtitle)
# fig.set_size_inches(7,7)
def on_resize(event):
    fig.tight_layout(rect=[0,0,1,1])
    fig.canvas.draw()

cid = fig.canvas.mpl_connect('resize_event', on_resize)
plt.show()
