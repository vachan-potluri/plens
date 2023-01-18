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
    "font.family": "Times",
    # "font.serif": ["Palatino"],
    "font.size": 10,
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



case_suffix = "chandrashekhar" # generally the flux scheme
res_dir = "result_16Jan2023_persson" # the name of result directory
sim_data_filename = "line_data.csv" # name of extracted line data file
# ref_soln_file = "../data/rho_ref_data.csv"
# ref_soln_file_delim = ","
ref_soln_file = "../data/rho_ref_data_monotonic.txt" # reference solution
ref_soln_file_delim = "\t"
plot_dir = "../plots/" # where to save the plots
N_values = [1,2,3,5]

ref_data = np.genfromtxt(ref_soln_file, delimiter=ref_soln_file_delim)
x_ref = ref_data[:,0]
rho_ref = ref_data[:,1]

# In case simulation data extracted has too much resolution, this will mask the resulting array
# to only include a fraction of the data
sim_sample_mask = np.arange(0, 6400, 4) # gives 1600 sampling points

# Zoomed ranges on x-axis: the x limits of the figure will be changed to obtain individual zoomed
# plots for the results. These figures will be stored as new figures
x_zooms_new = [ [0.3, 2.5], [-4,0] ]
x_zooms_new = []

# Range for the inset picture
x_range_inset = [0.3, 2.5]
# mask for plotting reference solution in inset
inset_mask_ref = np.logical_and(x_ref >= x_range_inset[0], x_ref <= x_range_inset[1])

for dof in [200,400,800]:
    fig, axes = plt.subplots(1,4,figsize=(10,2.5))
    for ax,N in zip(axes,N_values):
        sim_file = "dof{0}_{1}_{1}_N{2}_{3}/{4}/{5}".format(
            dof, N+1, N, case_suffix, res_dir, sim_data_filename
        )
        print("Reading {}".format(sim_file))
        sim_data = np.genfromtxt(sim_file, delimiter=",", names=True)
        x_sim = sim_data["Points0"][sim_sample_mask]
        rho_sim = sim_data["rho"][sim_sample_mask]
        ref_line, = ax.plot(
            x_ref,
            rho_ref,
            "b-",
            lw=1,
            label="Reference"
        )
        sim_line, = ax.plot(
            x_sim,
            rho_sim,
            "ro-",
            ms=0.5,
            lw=0.5,
            label="Simulation"
        )
        ax.grid()
        ax.set_title(r"$N={}$".format(N))

        # Now the inset
        # mask for plotting simulation data in inset
        inset_mask_sim = np.logical_and(x_sim >= x_range_inset[0], x_sim <= x_range_inset[1])
        ax_ins = ax.inset_axes([0.02,0.02,0.59,0.57]) # [0.02,0.02,0.59,0.57] works best
        ax_ins.plot(
            x_ref[inset_mask_ref],
            rho_ref[inset_mask_ref],
            "b-",
            lw=1,
        )
        ax_ins.plot(
            x_sim[inset_mask_sim],
            rho_sim[inset_mask_sim],
            "ro-",
            ms=0.5,
            lw=0.5
        )
        # remove all inset axes ticks and labels
        ax_ins.set_xticklabels([])
        ax_ins.set_yticklabels([])
        ax_ins.set_xticks([])
        ax_ins.set_yticks([])
        for s in ax_ins.spines.values(): s.set_edgecolor("darkgray") # inset border color
    # https://stackoverflow.com/a/6541454/8028836
    # plt.subplots_adjust(
    #     left=0.02,
    #     # bottom=None,
    #     right=0.99,
    #     top=0.79,
    #     wspace=0.1,
    #     # hspace=None
    # )
    # tight layout not working well for 1x4 array
    fig.tight_layout(
        pad=0.5
    ) # padding in fraction of font size
    # legend outside figure
    fig.legend(
        handles=[ref_line, sim_line],
        # labels=["Reference", "Simulation"],
        loc="lower center",
        bbox_to_anchor=(0.5,1),
        ncol=2, # number of columns
        borderpad=0.5,
        borderaxespad=0 # padding between axes and legend
    )
    plt.show()
    for fmt in ["png", "pdf"]:
        fig_filename = plot_dir + "dof{}_".format(dof) + case_suffix + ".{}".format(fmt)
        fig.savefig(fig_filename, format=fmt, bbox_inches="tight", pad_inches=2./72)
        print("Saved figure {}".format(fig_filename))
    
    for i,x_lim in enumerate(x_zooms_new):
        for row in axes:
            for ax in row:
                ax.set_xlim(x_lim)
                autoscale_y(ax)
        # fig.set_size_inches(7.5, 3)
        fig.tight_layout() # don't pass any arguments for second time
        # plt.show() # doesn't show second time
        for fmt in ["png", "pdf"]:
            fig_filename = (plot_dir + "dof{}_".format(dof) + case_suffix +
                "_zoom{}.{}".format(i+1,fmt))
            fig.savefig(fig_filename, format=fmt)
            print("Saved figure {}".format(fig_filename))
