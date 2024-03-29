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

parser = argparse.ArgumentParser(
    description = "A script to plot extracted data for Kitamura's case."
)
parser.add_argument("res_dir", help="Result directory (absolute or relative).")
parser.add_argument("counter", help="Output counter.", type=int)
parser.add_argument(
    "-a",
    "--alpha",
    help="Plot alpha in heat transfer figure",
    action="store_true"
)
args = parser.parse_args()

res_dir = args.res_dir
if res_dir[-1] != "/":
    res_dir += "/"
counter = args.counter

surface_data = np.genfromtxt(
    "{}surface_data_{}.csv".format(res_dir, counter), delimiter=",", names=True
)
kitamura_p_data = np.genfromtxt("../data/kitamura_p_ausmpw_coarse.csv", delimiter=",")
kitamura_q_ausmpw_data = np.genfromtxt("../data/kitamura_q_ausmpw_coarse.csv", delimiter=",")
kitamura_q_roe_data = np.genfromtxt("../data/kitamura_q_roe_efix_coarse.csv", delimiter=",")

r = 0.02
p_inf = 370.7
p0 = p_inf*84.9 # this ratio is given in Kitamura (2010)
q_fr = 17.5e4 # given in Kitamura et al (2010): Fay & Riddells value
p = surface_data["p"]
theta = np.arcsin(surface_data["Points1"]/r)
theta_deg = theta*180/np.pi
alpha = surface_data["alpha"]
q = surface_data["qx"]*np.cos(theta) - surface_data["qy"]*np.sin(theta)

fig, ax = plt.subplots(1,1)
ax.plot(theta_deg, p/p0, "r-", label="PLENS")
ax.plot(kitamura_p_data[:,0], kitamura_p_data[:,1], "b--", label="Kitamura et al\n(2010) simulation")
ax.legend()
ax.grid()
ax.set_xlabel(r"$\theta$ [degrees]")
ax.set_ylabel(r"$\displaystyle\frac{p}{p_{10}}$", rotation=0, labelpad=20)
fig.tight_layout(pad=0.25)
plt.show()
for fmt in ["png", "pdf"]:
    fig_filename = "{}p_{}.{}".format(res_dir, counter, fmt)
    fig.savefig(fig_filename, format=fmt)
    print("Saved figure {}".format(fig_filename))
del fig, ax

fig, ax = plt.subplots(1,1)
# mask = np.abs(theta_deg) < 50
# ax.plot(theta_deg[mask], q[mask]/q_fr, "r-", label="PLENS")
ax.plot(theta_deg, q/q_fr, "r-", label="PLENS")
ax.plot(kitamura_q_ausmpw_data[:,0], kitamura_q_ausmpw_data[:,1], "b--", label="Kitamura et al (2010)\nsimulation (AUSMPW+)")
ax.plot(kitamura_q_roe_data[:,0], kitamura_q_roe_data[:,1], "g-.", label="Kitamura et al (2010)\nsimulation (Roe E-fix)")
if args.alpha: ax.plot(theta*180/np.pi, alpha, c="magenta", ls=":", label=r"$\alpha$")
ax.legend()
ax.grid()
ax.set_xlabel(r"$\theta$ [degrees]")
ax.set_ylabel(r"$\displaystyle\frac{q^{''}}{q^{''}_{\textrm{FR}}}$", rotation=0, labelpad=20)
fig.tight_layout(pad=0.25)
plt.show()
for fmt in ["png", "pdf"]:
    fig_filename = "{}q_{}.{}".format(res_dir, counter, fmt)
    fig.savefig(fig_filename, format=fmt)
    print("Saved figure {}".format(fig_filename))
