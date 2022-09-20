# Plots group plots using saved data from plot_individual.py for a given dof

import numpy as np
import matplotlib.pyplot as plt
import argparse
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 10,
    "axes.formatter.limits": [-2,2],
})

parser = argparse.ArgumentParser(
    description = "A script to compare numerical and exact results for different values of N, for "
    "a given number of dof and flux scheme. It uses the data saved by 'plot_individual.py' for "
    "this purpose. The location data files will be assumed to be of the following form:\n\n"
    "\t[case_dir]/dof[dof]_12_12_N1_[flux]/[result_dir]/[data_filename]"
    "\n\nfor N=1 and so on for N=2,3,5."
)
parser.add_argument(
    "case_dir",
    help="The case directory. This must be the 'run' directory from which individual simulation "
    "folders can be directly accessed."
)
parser.add_argument(
    "dof",
    type=int,
    help="Degrees of freedom in x-direction. Presently y and z directional dofs are taken as 12 "
    "each."
)
parser.add_argument(
    "flux",
    help="The flux scheme used."
)
parser.add_argument(
    "result_dir",
    help="The result directory."
)
parser.add_argument(
    "data_filename",
    help="The filename of individual data files written by 'plot_individual.py'."
)
parser.add_argument(
    "plot_title",
    help="The title to be given to the group plot. Latex expression can be given."
)
parser.add_argument(
    "-s",
    "--save",
    help="This option can be used to save the generated plot. The first argument will be the "
    "directory, and the second argument will be the base file name. The plot will be saved in "
    "'pdf' and 'png' formats.",
    nargs=2,
    default=["", ""],
    action="store"
)
parser.add_argument(
    "-z",
    "--size",
    help="The size of the plot window (in inches) to be used.",
    nargs=2,
    type=float,
    default=[7,7],
    action="store"
)
args = parser.parse_args()

N_values = [1,2,3,5]
subtitles = []
data_files = []
for N in N_values:
    subtitles.append(r"$N={}$".format(N))
    data_files.append("{}/dof{}_12_12_N{}_{}/{}/{}".format(
        args.case_dir,
        args.dof,
        N,
        args.flux,
        args.result_dir,
        args.data_filename
    ))
figtitle = r"{}".format(args.plot_title)



# Override data_files
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
# data_files = [
#     "../test1-4/run/dof400_12_12_N1_hllc/result/comparison_data.csv",
#     "../test1-4/run/dof400_12_12_N2_hllc/result/comparison_data.csv",
#     "../test1-4/run/dof400_12_12_N3_hllc/result/comparison_data.csv",
#     "../test1-4/run/dof400_12_12_N5_hllc/result/comparison_data.csv"
# ]
# data_files = [
#     "../test1-5/run/dof400_12_12_N1_rusanov/result/comparison_data.csv",
#     "../test1-5/run/dof400_12_12_N2_rusanov/result/comparison_data.csv",
#     "../test1-5/run/dof400_12_12_N3_rusanov/result/comparison_data.csv",
#     "../test1-5/run/dof400_12_12_N5_rusanov/result/comparison_data.csv"
# ]
# data_files = [
#     "../test1-1/run/dof800_12_12_N1_hllc/result/comparison_data.csv",
#     "../test1-1/run/dof800_12_12_N2_hllc/result/comparison_data.csv",
#     "../test1-1/run/dof800_12_12_N3_hllc/result/comparison_data.csv",
#     "../test1-1/run/dof800_12_12_N5_hllc/result1/comparison_data.csv"
# ]
# data_files = [
#     "../test1-2/run/dof800_12_12_N1_rusanov/result/comparison_data.csv",
#     "../test1-2/run/dof800_12_12_N2_rusanov/result/comparison_data.csv",
#     "../test1-2/run/dof800_12_12_N3_rusanov/result/comparison_data.csv",
#     "../test1-2/run/dof800_12_12_N5_rusanov/result/comparison_data.csv"
# ]
# data_files = [
#     "../test1-3/run/dof800_12_12_N1_hllc/result/comparison_data.csv",
#     "../test1-3/run/dof800_12_12_N2_hllc/result/comparison_data.csv",
#     "../test1-3/run/dof800_12_12_N3_hllc/result/comparison_data.csv",
#     "../test1-3/run/dof800_12_12_N5_hllc/result1/comparison_data.csv"
# ]
# data_files = [
#     "../test1-4/run/dof800_12_12_N1_hllc/result/comparison_data.csv",
#     "../test1-4/run/dof800_12_12_N2_hllc/result/comparison_data.csv",
#     "../test1-4/run/dof800_12_12_N3_hllc/result/comparison_data.csv",
#     "../test1-4/run/dof800_12_12_N5_hllc/result1/comparison_data.csv"
# ]
# data_files = [
#     "../test1-5/run/dof800_12_12_N1_rusanov/result/comparison_data.csv",
#     "../test1-5/run/dof800_12_12_N2_rusanov/result/comparison_data.csv",
#     "../test1-5/run/dof800_12_12_N3_rusanov/result/comparison_data.csv",
#     "../test1-5/run/dof800_12_12_N5_rusanov/result1/comparison_data.csv"
# ]

fig, axes = plt.subplots(2,2)
i = 0
for row in axes:
    for ax in row:
        print("Reading data from {}".format(data_files[i]))
        data = np.genfromtxt(data_files[i], delimiter=",")
        ax.plot(data[:,0], data[:,1], "b--", lw=1, label="Exact")
        ax.plot(data[:,0], data[:,2], "r-", lw=1, label="Simulation")
        ax.set_title(subtitles[i])
        ax.grid()
        if i==0: ax.legend(loc="best", handlelength=1)
        i += 1
if figtitle != "": fig.suptitle(figtitle)
fig.set_size_inches(args.size[0], args.size[1])
fig.tight_layout(rect=[0,0,1,1], pad=0.75)
# def on_resize(event):
#     fig.tight_layout(rect=[0,0,1,1])
#     fig.canvas.draw()
# cid = fig.canvas.mpl_connect('resize_event', on_resize)
plt.show()

if args.save[0] != "" and args.save[1] != "":
    for fmt in ["png", "pdf"]:
        full_plot_filename = "{}/{}.{}".format(args.save[0], args.save[1], fmt)
        fig.savefig(full_plot_filename, format=fmt)
        print("Written file {}".format(full_plot_filename))
