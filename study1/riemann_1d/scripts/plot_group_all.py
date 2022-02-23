# Plots group plots using saved data from plot_individual.py for all dofs

import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 10
plt.rcParams["mathtext.fontset"] = "dejavuserif"

parser = argparse.ArgumentParser(
    description = "A script to compare numerical and exact results for different values of N, for "
    "all dofs with a given flux scheme. It uses the data saved by 'plot_individual.py' for "
    "this purpose. The location data files will be assumed to be of the following form:\n\n"
    "\t[case_dir]/dof200_12_12_N1_[flux]/[result_dir]/[data_filename]"
    "\n\nfor 200 dofs and N=1, and so on for 400, 800 dofs and N=2,3,5."
)
parser.add_argument(
    "case_dir",
    help="The case directory. This must be the 'run' directory from which individual simulation "
    "folders can be directly accessed."
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
    type=int,
    default=[7,7],
    action="store"
)
args = parser.parse_args()

N_values = [1,2,3,5]
dofs = [200, 400, 800]
subtitles = pd.Series(data="", index=N_values)
data_files = pd.DataFrame(data="", index=N_values, columns=dofs)
for N in N_values:
    subtitles.loc[N] = r"$N={}$".format(N)
    for dof in dofs:
        data_files.loc[N, dof] = "{}/dof{}_12_12_N{}_{}/{}/{}".format(
            args.case_dir,
            dof,
            N,
            args.flux,
            args.result_dir,
            args.data_filename
        )
figtitle = r"{}".format(args.plot_title)
dof_linestyles = pd.Series(["-", "--", "-."], index=dofs)
dof_colors = pd.Series(["r", "g", "b"], index=dofs)

fig, axes = plt.subplots(2,2)
i = 0
for row in axes:
    for ax in row:
        N = N_values[i]
        for dof in dofs:
            print("Reading data from {}".format(data_files.loc[N,dof]))
            data = np.genfromtxt(data_files.loc[N,dof], delimiter=",")
            # plot exact data in dof 200 case
            if dof == 200:
                ax.plot(data[:,0], data[:,1], "k:", label="Exact")
            ax.plot(
                data[:,0],
                data[:,2],
                ls=dof_linestyles[dof],
                c=dof_colors[dof],
                alpha=0.75,
                label="{} DoFs".format(dof)
            )
        ax.set_title(subtitles[N])
        ax.grid()
        ax.legend(loc="best")
        i += 1
fig.suptitle(figtitle)
fig.set_size_inches(args.size[0], args.size[1])
fig.tight_layout(rect=[0,0,1,1])
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
