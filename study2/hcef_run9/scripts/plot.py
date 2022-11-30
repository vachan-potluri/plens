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

p_inf = 30.3845
rho_inf = 0.000845
u_inf = 2566.0
L = 10.17e-2

parser = argparse.ArgumentParser(
    description = "A script to plot extracted data for HCEF case."
)
parser.add_argument("res_dir", help="Result directory (absolute or relative).")
parser.add_argument("counter", help="Output counter.", type=int)
parser.add_argument(
    "-a",
    "--azimuths",
    help="Extract data on multiple azimuthal angles (in degrees). Default: ''. If this argument is "
        + "empty, only the zero-degree azimuthal angle is used for extracting.",
    default=[0.0],
    type=float,
    action="store",
    nargs="*"
)
args = parser.parse_args()

res_dir = args.res_dir
if res_dir[-1] != "/":
    res_dir += "/"
counter = args.counter
harvey_cp = np.genfromtxt("../data/cp_harvey.csv", delimiter=",")
gnoffo_cp = np.genfromtxt("../data/cp_gnoffo.csv", delimiter=",")
harvey_St = np.genfromtxt("../data/St_harvey.csv", delimiter=",")
gnoffo_St = np.genfromtxt("../data/St_gnoffo.csv", delimiter=",")

fig1, ax1 = plt.subplots(1,1)
fig2, ax2 = plt.subplots(1,1)
for azimuth in args.azimuths:
    phi = np.pi*azimuth/180
    cylinder_data = np.genfromtxt(
        "{}cylinder_data_phi{:.0f}_{}.csv".format(res_dir, azimuth, counter),
        delimiter=",",
        names=True
    )
    flare_data = np.genfromtxt(
        "{}flare_data_phi{:.0f}_{}.csv".format(res_dir, azimuth, counter),
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
    ax1.plot(x[cp_mask], cp[cp_mask], ls="--", label=r"$\phi={:.0f}$".format(azimuth), alpha=0.5)



    St_cylinder = -2*cylinder_data["qy"]/(rho_inf*u_inf**3)
    St_flare = -2*(
        (np.cos(phi)*flare_data["qy"] + np.sin(phi)*flare_data["qz"])*np.cos(np.pi/6) -
        flare_data["qx"]*np.sin(np.pi/6)
    )/(rho_inf*u_inf**3)
    St = np.concatenate([St_cylinder, St_flare])
    St_mask = np.isfinite(St)
    ax2.plot(x[St_mask], St[St_mask], ls="--", label=r"$\phi={:.0f}$".format(azimuth), alpha=0.5)

ax1.plot(harvey_cp[:,0], harvey_cp[:,1], "bo", label="Harvey et al (2001)\nexperiment")
ax1.plot(gnoffo_cp[:,0], gnoffo_cp[:,1], "g^", label="Gnoffo (2001)\nsimulation")
ax1.set_xlabel(r"$x/L$")
ax1.set_ylabel(r"$\displaystyle\frac{p-p_\infty}{\frac{1}{2}\rho_\infty u_\infty^2}$", rotation=0, labelpad=25)
ax1.grid()
ax1.legend(loc="best")
fig1.tight_layout(pad=0.25)
for fmt in ["png", "pdf"]:
    full_figname = "{}cp_comparison_{}.{}".format(res_dir, counter, fmt)
    fig1.savefig(full_figname, format=fmt)
    print("Written file {}".format(full_figname))

ax2.plot(harvey_St[:,0], harvey_St[:,1], "bo", label="Harvey et al (2001)\nexperiment")
ax2.plot(gnoffo_St[:,0], gnoffo_St[:,1], "g^", label="Gnoffo (2001)\nsimulation")
ax2.set_xlabel(r"$x/L$")
ax2.set_ylabel(r"$\displaystyle\frac{2q^{\prime\prime}_w}{\rho_\infty u_\infty^3}$", rotation=0, labelpad=20)
ax2.grid()
ax2.legend(loc="best")
fig2.tight_layout(pad=0.25)
for fmt in ["png", "pdf"]:
    full_figname = "{}St_comparison_{}.{}".format(res_dir, counter, fmt)
    fig2.savefig(full_figname, format=fmt)
    print("Written file {}".format(full_figname))
plt.show()



"""
#del fig, ax
fig, ax = plt.subplots(1,1)
ax.plot(x[txy_mask], txy[txy_mask], "r-", label="PLENS")
# Gnoffo's results on unadaptive 272 by 96 grid
x_sep_gnoffo = 0.52
x_att_gnoffo = 1.31
txy_max = np.max(txy[txy_mask])
txy_min = np.min(txy[txy_mask])
txy_range = np.array([txy_min, -txy_min])
ax.plot(x_sep_gnoffo*np.ones(2), txy_range, "g--")
ax.plot(x_att_gnoffo*np.ones(2), txy_range, "g--", label="Gnoffo (2001) separation\n& attachment locations")
ax.legend()
ax.minorticks_on()
ax.grid(which="major", linestyle='-', linewidth='0.5')
ax.grid(which="minor", linestyle='-', linewidth='0.25')
ax.set_xlabel(r"$x/L$")
ax.set_ylabel(r"$\tau_{xy}$ [Pa]")
fig.tight_layout()
fig.savefig("{}txy_comparison_{}.png".format(res_dir, counter), format="png")
fig.savefig("{}txy_comparison_{}.pdf".format(res_dir, counter), format="pdf")
plt.show()
"""
