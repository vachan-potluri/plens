import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    # "font.serif": ["Palatino"],
    "font.size": 14,
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



ref_data = np.genfromtxt("../data/rho_ref_data_monotonic.txt", delimiter="\t")
sim_data = np.genfromtxt(
    "../run/dof200_2_2_N1_chandrashekhar/result_16Jun2022/line_data.csv",
    delimiter=",",
    names=True
)
sim_sample_mask = np.arange(0, 6400, 4)

fig, ax = plt.subplots(1,1)
ax.plot(
    ref_data[:,0],
    ref_data[:,1],
    "b-",
    lw=1,
    label="Reference"
)
ax.plot(
    sim_data["Points0"][sim_sample_mask],
    sim_data["rho"][sim_sample_mask],
    "ro-",
    ms=1,
    lw=0.5,
    label="Simulation"
)
ins_ref_mask = np.logical_and(ref_data[:,0] > 0.3, ref_data[:,0] < 2.5)
ins = ax.inset_axes([0.01,0.01,0.58,0.58])
ins.plot(ref_data[:,0][ins_ref_mask], ref_data[:,1][ins_ref_mask])
ins.set_xticklabels([])
ins.set_yticklabels([])
ins.set_xticks([])
ins.set_yticks([])
for s in ins.spines.values(): s.set_edgecolor("darkgray")
ax.legend()
plt.show()
