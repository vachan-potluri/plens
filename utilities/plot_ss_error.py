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
mask = (ss_error != 1) # indicates simulation start and/or restart
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
plt.show()
for fmt in ["png", "pdf"]:
    figname = "ss_error_plot.{}".format(fmt)
    fig.savefig(figname)
    print("Saved figure {}".format(figname))
