import argparse
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
    "font.size": 12,
    "axes.formatter.limits": [-2,2]
})
import pandas as pd

p_inf = 30.3845
rho_inf = 0.000845
u_inf = 2566.0
L = 10.17e-2

ids = ["Chandrashekhar", "Rusanov-HLLC"]
res_dirs = pd.Series(
    data=[
        "../trial2/result_15Jun2022/",
        "../trial2/result2_15Jun2022/"
    ],
    index=ids
)
res_ctrs = pd.Series(data=[1873, 1938], index=ids)
labels = pd.Series(data=["Chandrashekhar", "Rusanov-HLLC"], index=ids)
plot_formats = pd.Series(data=["r-", "m--"], index=ids)
plot_dir = "../plots/trial2/"
fig_basename = "result_15Jun2022_result2_15Jun2022"



harvey_cp = np.genfromtxt("../data/cp_harvey.csv", delimiter=",")
gnoffo_cp = np.genfromtxt("../data/cp_gnoffo.csv", delimiter=",")
harvey_St = np.genfromtxt("../data/St_harvey.csv", delimiter=",")
gnoffo_St = np.genfromtxt("../data/St_gnoffo.csv", delimiter=",")

fig1, ax1 = plt.subplots(1,1)
fig2, ax2 = plt.subplots(1,1)

ax1.plot(harvey_cp[:,0], harvey_cp[:,1], "bo", label="Harvey et al (2001)\nexperiment")
ax1.plot(gnoffo_cp[:,0], gnoffo_cp[:,1], "g^", label="Gnoffo (2001)\nsimulation")
ax2.plot(harvey_St[:,0], harvey_St[:,1], "bo", label="Harvey et al (2001)\nexperiment")
ax2.plot(gnoffo_St[:,0], gnoffo_St[:,1], "g^", label="Gnoffo (2001)\nsimulation")

for idx in ids:
    cylinder_data = np.genfromtxt(
        "{}/cylinder_data_{}.csv".format(res_dirs.loc[idx], res_ctrs.loc[idx]),
        delimiter=",",
        names=True
    )
    flare_data = np.genfromtxt(
        "{}/flare_data_{}.csv".format(res_dirs.loc[idx], res_ctrs.loc[idx]),
        delimiter=",",
        names=True
    )
    x_cylinder = cylinder_data["Points0"]/L
    x_flare = flare_data["Points0"]/L
    x = np.concatenate([x_cylinder, x_flare])

    cp_cylinder = 2*(cylinder_data["p"] - p_inf)/(rho_inf*u_inf**2)
    cp_flare = 2*(flare_data["p"] - p_inf)/(rho_inf*u_inf**2)
    cp = np.concatenate([cp_cylinder, cp_flare])
    cp_mask = np.isfinite(cp)

    St_cylinder = -2*cylinder_data["qy"]/(rho_inf*u_inf**3)
    St_flare = -2*(flare_data["qy"]*np.cos(np.pi/6) - flare_data["qx"]*np.sin(np.pi/6))/(rho_inf*u_inf**3)
    St = np.concatenate([St_cylinder, St_flare])
    St_mask = np.isfinite(St)

    ax1.plot(x[cp_mask], cp[cp_mask], plot_formats.loc[idx], label=labels.loc[idx])
    ax2.plot(x[St_mask], St[St_mask], plot_formats.loc[idx], label=labels.loc[idx])

ax1.set_xlabel(r"$x/L$")
ax1.set_ylabel(r"$\displaystyle\frac{p-p_\infty}{\frac{1}{2}\rho_\infty u_\infty^2}$", rotation=0, labelpad=25)
ax1.grid()
ax1.legend(loc="best")
fig1.tight_layout(pad=0.25)
ax2.set_xlabel(r"$x/L$")
ax2.set_ylabel(r"$\displaystyle\frac{2q^{\prime\prime}_w}{\rho_\infty u_\infty^3}$", rotation=0, labelpad=20)
ax2.grid()
ax2.legend(loc="best")
fig2.tight_layout(pad=0.25)

plt.show()

for fmt in ["png", "pdf"]:
    fig1_name = "{}/{}_cp.{}".format(plot_dir, fig_basename, fmt)
    fig2_name = "{}/{}_St.{}".format(plot_dir, fig_basename, fmt)
    fig1.savefig(fig1_name, format=fmt)
    fig2.savefig(fig2_name, format=fmt)
    print("Written {} and {}".format(fig1_name, fig2_name))
