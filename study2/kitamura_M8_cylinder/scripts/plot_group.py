import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
    "font.size": 14,
    "axes.formatter.limits": [-2,2]
})

kitamura_p_data = np.genfromtxt("../data/kitamura_p_ausmpw_coarse.csv", delimiter=",")
kitamura_q_ausmpw_data = np.genfromtxt("../data/kitamura_q_ausmpw_coarse.csv", delimiter=",")
kitamura_q_roe_data = np.genfromtxt("../data/kitamura_q_roe_efix_coarse.csv", delimiter=",")

r = 0.02
p_inf = 370.7
p0 = p_inf*84.9 # this ratio is given in Kitamura (2010)
q_fr = 17.5e4 # given in Kitamura et al (2010): Fay & Riddells value

ids = ["chandrashekhar2", "rusanov-hllc2", "chandrashekhar1"]
res_dirs = pd.Series(
    ["../trial6/restart_22Jun2022", "../trial6/result2_01Jul2022", "../trial7/result_30Jun2022"],
    index=ids
)
counters = pd.Series([1206, 650, 437], index=ids)
labels = pd.Series(
    [r"$N=2$ Chandrashekhar", r"$N=2$ Rusanov-HLLC", r"$N=1$ Chandrashekhar"],
    index=ids
)
plot_styles = pd.Series(["r-", "b--", "g-."], index=ids)

fig, ax = plt.subplots(1,1)
theta_deg_ref = kitamura_q_ausmpw_data[:,0]
q_ref = kitamura_q_ausmpw_data[:,1]
mask_ref = np.abs(theta_deg_ref) < 40
ax.plot(
    theta_deg_ref[mask_ref],
    q_ref[mask_ref],
    "m:",
    label="Kitamura et al (2010)\nsimulation (AUSMPW+)"
)
for idx in ids:
    surface_data = np.genfromtxt(
        "{}/surface_data_{}.csv".format(res_dirs.loc[idx], counters.loc[idx]),
        delimiter=",",
        names=True
    )
    theta = np.arcsin(surface_data["Points1"]/r)
    theta_deg = theta*180/np.pi
    mask = np.abs(theta_deg) <= 40
    q = surface_data["qx"]*np.cos(theta) - surface_data["qy"]*np.sin(theta)
    ax.plot(theta_deg[mask], q[mask]/q_fr, plot_styles.loc[idx], label=labels.loc[idx])
ax.legend()
ax.grid()
ax.set_xlabel(r"$\theta$ [degrees]")
ax.set_ylabel(
    r"$\displaystyle\frac{q^{\prime\prime}}{q^{\prime\prime}_{\textrm{FR}}}$",
    rotation=0,
    labelpad=20
)
fig.tight_layout(pad=0.25)
plt.show()
for fmt in ["png", "pdf"]:
    fig_filename = "../plots/q_group.{}".format(fmt)
    fig.savefig(fig_filename, format=fmt)
    print("Saved figure {}".format(fig_filename))
