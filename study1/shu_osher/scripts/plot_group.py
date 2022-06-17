# Script for group plotting
# Introduced on 17-Jun-2022 for SIS2022 analysis

# This must be executed from the "run" directory from where individual case directories are
# accessible. The result line data will be assumed to have the same name in all directories. The
# plots will be saved in "../plots" directory (note that this script will be executed from run
# directory).

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    # "font.serif": ["Palatino"],
    "font.size": 12,
    # "axes.formatter.limits": [-2,2]
})



# https://stackoverflow.com/a/35094823/8028836
def autoscale_y(ax,margin=0.1):
    """This function rescales the y-axis based on the data that is visible given the current xlim of the axis.
    ax -- a matplotlib axes object
    margin -- the fraction of the total height of the y-data to pad the upper and lower ylims"""

    import numpy as np

    def get_bottom_top(line):
        xd = line.get_xdata()
        yd = line.get_ydata()
        lo,hi = ax.get_xlim()
        y_displayed = yd[((xd>lo) & (xd<hi))]
        h = np.max(y_displayed) - np.min(y_displayed)
        bot = np.min(y_displayed)-margin*h
        top = np.max(y_displayed)+margin*h
        return bot,top

    lines = ax.get_lines()
    bot,top = np.inf, -np.inf

    for line in lines:
        new_bot, new_top = get_bottom_top(line)
        if new_bot < bot: bot = new_bot
        if new_top > top: top = new_top

    ax.set_ylim(bot,top)



case_suffix = "chandrashekhar"
res_dir = "result_16Jun2022"
sim_data_filename = "line_data.csv"
# ref_soln_file = "../data/rho_ref_data.csv"
# ref_soln_file_delim = ","
ref_soln_file = "../data/rho_ref_data_monotonic.txt"
ref_soln_file_delim = "\t"
plot_dir = "../plots/"
N_values = [1,2,3,5]

ref_data = np.genfromtxt(ref_soln_file, delimiter=ref_soln_file_delim)

# In case simulation data extracted has too much resolution, this will mask the resulting array
# to only include a fraction of the data
sim_sample_mask = np.arange(0, 6400, 4) # gives 1600 sampling points

# Zoomed ranges on x-axis: the x limits of the figure will be changed to obtain individual zoomed
# plots for the results
x_zooms = [ [0.3, 2.5], [-3,-1] ]

for dof in [200,400,800]:
    fig, axes = plt.subplots(2,2)
    N_id = 0
    for row in axes:
        for ax in row:
            N = N_values[N_id]
            sim_file = "dof{0}_{1}_{1}_N{2}_{3}/{4}/{5}".format(
                dof, N+1, N, case_suffix, res_dir, sim_data_filename
            )
            print("Reading {}".format(sim_file))
            sim_data = np.genfromtxt(sim_file, delimiter=",", names=True)
            ax.plot(
                ref_data[:,0],
                ref_data[:,1],
                "b-",
                lw=1,
                label="Reference")
            ax.plot(
                sim_data["Points0"][sim_sample_mask],
                sim_data["rho"][sim_sample_mask],
                "ro-",
                ms=1,
                lw=1,
                label="Simulation")
            ax.grid()
            ax.legend(handlelength=1)
            ax.set_title(r"$N={}$".format(N))
            N_id += 1
    fig.tight_layout(rect=[0,0,1,1], pad=0.25)
    fig.set_size_inches(7.5, 6)
    plt.show()
    for fmt in ["png", "pdf"]:
        fig_filename = plot_dir + "dof{}_".format(dof) + case_suffix + ".{}".format(fmt)
        fig.savefig(fig_filename, format=fmt)
        print("Saved figure {}".format(fig_filename))
    
    for i,x_lim in enumerate(x_zooms):
        for row in axes:
            for ax in row:
                ax.set_xlim(x_lim)
                autoscale_y(ax)
        fig.tight_layout() # don't pass any arguments for second time
        # plt.show() # doesn't show second time
        for fmt in ["png", "pdf"]:
            fig_filename = (plot_dir + "dof{}_".format(dof) + case_suffix +
                "_zoom{}.{}".format(i+1,fmt))
            fig.savefig(fig_filename, format=fmt)
            print("Saved figure {}".format(fig_filename))
