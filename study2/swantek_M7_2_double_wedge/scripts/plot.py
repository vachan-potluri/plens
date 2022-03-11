import argparse
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams["font.size"] = 10
plt.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}"
# plt.rcParams["figure.figsize"] = 7,6

l_wedge1 = 50.8e-3
l_wedge2 = 25.4e-3
theta_wedge1 = 30*np.pi/180
theta_wedge2 = 55*np.pi/180

parser = argparse.ArgumentParser(
    description = "A script to plot extracted data for Swantek's double wedge geometry"
)
parser.add_argument("res_dir", help="Result directory (absolute or relative).")
parser.add_argument("counter", help="Output counter.", type=int)
args = parser.parse_args()

res_dir = args.res_dir
if res_dir[-1] != "/":
    res_dir += "/"
counter = args.counter
wedge1_data = np.genfromtxt(
    "{}wedge1_data_{}.csv".format(res_dir, counter), delimiter=",", names=True
)
wedge2_data = np.genfromtxt(
    "{}wedge2_data_{}.csv".format(res_dir, counter), delimiter=",", names=True
)
komives_q = np.genfromtxt("../data/komives_q_0.27ms.csv", delimiter=",")

x_cylinder = wedge1_data["Points0"]/l_wedge1
x_flare = wedge2_data["Points0"]/l_wedge1
x = np.concatenate([x_cylinder, x_flare])

q_wedge1 = wedge1_data["qx"]*np.sin(theta_wedge1) - wedge1_data["qy"]*np.cos(theta_wedge1)
q_wedge2 = wedge2_data["qx"]*np.sin(theta_wedge2) - wedge2_data["qy"]*np.cos(theta_wedge2)
q = np.concatenate([q_wedge1, q_wedge2])
q_mask = np.isfinite(q)

fig, ax = plt.subplots(1,1)
ax.plot(
    komives_q[:,0], komives_q[:,1], "bo", markersize=4,
    label="Komives et al (2014)\nsimulation at \SI{0.27}{\milli\second}"
)
ax.plot(x[q_mask], q[q_mask], "r-", label="PLENS")
ax.set_xlabel(r"$x/L$")
ax.set_ylabel(r"$q^{\prime\prime}_w$")
ax.grid()
ax.legend(loc="best")
fig.tight_layout()
for fmt in ["png", "pdf"]:
    full_figname = "{}heat_flux_comparison_{}.{}".format(res_dir, counter, fmt)
    fig.savefig(full_figname, format=fmt)
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
