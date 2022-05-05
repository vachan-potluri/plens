# A small script to plot the error data from `./output.ss_error`
# The file must be available in the working directory

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
    "font.size": 12,
    "axes.formatter.limits": [-2,2]
})

data = np.genfromtxt("output.ss_error")
counter = data[:,0]
time = data[:,1]
ss_error = data[:,2]
mask = np.ones(ss_error.size, dtype=bool) # indicates simulation start and/or restart
for i in range(mask.size):
    if i == 0: mask[i] = False # first value will always be masked, since it is a start of simulation
    if i != mask.size-1:
        if counter[i] == counter[i+1]:
            # indicates restart, mask the restarted values
            mask[i+1] = False
fig, axes = plt.subplots(2,1)
ax = axes[0]
ax.plot(time[mask], ss_error[mask], "b-")
ax.set_yscale("log")
ax.set_xlabel(r"$t$ [sec]")
ax.set_ylabel("Residual")
ax.grid()
ax = axes[1]
ax.plot(counter[mask], ss_error[mask], "b-")
ax.set_yscale("log")
ax.set_xlabel("Output counter")
ax.set_ylabel("Residual")
ax.grid()
fig.tight_layout()
fig.set_size_inches(6,9)
plt.show()
for fmt in ["png", "pdf"]:
    figname = "ss_error_plot.{}".format(fmt)
    fig.savefig(figname)
    print("Saved figure {}".format(figname))
