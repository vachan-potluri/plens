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

## Some important notes
# cp is p/rho u^2
# Although the paper says cp is taken as p/(0.5 rho u^2), the values don't match when this
# expression is used. I have checked the freestream pressure and verified that the cp formula
# is missing a factor 2

# Also, the paper plots quantities vs wall coordinate, instead of x-coordinate

p_inf = 9.99473
rho_inf = 0.0004827
u_inf = 2401.824
T_inf = 72.22
T_w = 297.22
cp_air = 1005
L = 0.438912
wedge_angle = 24*np.pi/180

parser = argparse.ArgumentParser(
    description = "A script to plot extracted data for Hung wedge 18 case."
)
parser.add_argument("res_dir", help="Result directory (absolute or relative).")
parser.add_argument("counter", help="Output counter.", type=int)
args = parser.parse_args()

res_dir = args.res_dir
if res_dir[-1] != "/":
    res_dir += "/"
counter = args.counter
plate_data = np.genfromtxt(
    "{}plate_data_{}.csv".format(res_dir, counter), delimiter=",", names=True
)
ramp_data = np.genfromtxt(
    "{}ramp_data_{}.csv".format(res_dir, counter), delimiter=",", names=True
)

holden_cp = np.genfromtxt("../data/holden_cp.csv", delimiter=",")
hung_cp = np.genfromtxt("../data/hung_cp.csv", delimiter=",")
rudy_cp_2d = np.genfromtxt("../data/rudy_cp_2d_201.csv", delimiter=",") # log10(50 cp), x/L
# holden_ch = np.genfromtxt("../data/holden_ch.csv", delimiter=",")
# hung_ch = np.genfromtxt("../data/hung_ch.csv", delimiter=",")
# holden_cf = np.genfromtxt("../data/holden_cf.csv", delimiter=",")
# hung_cf = np.genfromtxt("../data/hung_cf.csv", delimiter=",")

x_plate = plate_data["Points0"]/L
x_ramp = ramp_data["Points0"]/L
s_wall = np.concatenate([x_plate, 1+(x_ramp-1)/np.cos(wedge_angle)])
x_wall = np.concatenate([x_plate, x_ramp])

cp_plate = 2*(plate_data["p"])/(rho_inf*u_inf**2)
cp_ramp = 2*(ramp_data["p"])/(rho_inf*u_inf**2)
cp = np.concatenate([cp_plate, cp_ramp])
cp_mask = np.isfinite(cp)

# Holden's data and Hung's data both use wall distance, transforming it to x-coordinate
x_holden = holden_cp[:,0]
for i in range(x_holden.size):
    if x_holden[i] > 1:
        x_holden[i] = 1 + (x_holden[i]-1)*np.cos(wedge_angle)
x_hung = hung_cp[:,0]
for i in range(x_holden.size):
    if x_hung[i] > 1:
        x_hung[i] = 1 + (x_hung[i]-1)*np.cos(wedge_angle)

# ch_plate = -plate_data["qy"]/(rho_inf*u_inf*(0.5*u_inf**2 + cp_air*(T_inf-T_w)))
# ch_ramp = -(
#     ramp_data["qy"]/np.cos(wedge_angle) -
#     ramp_data["qx"]*np.sin(wedge_angle)
# )/(rho_inf*u_inf*(0.5*u_inf**2 + cp_air*(T_inf-T_w)))
# ch = np.concatenate([ch_plate, ch_ramp])
# ch_mask = np.isfinite(ch)

# cf_plate = 2*plate_data["txy"]/(rho_inf*u_inf**2)
# cf_ramp = 2*(
#     ramp_data["txy"]
#     # ramp_data["txy"]*np.cos(wedge_angle) + ramp_data["tyy"]*np.sin(wedge_angle) +
#     # ramp_data["txx"]*np.cos(wedge_angle) + ramp_data["txy"]*np.sin(wedge_angle)
#     )/(rho_inf*u_inf**2)
# cf = np.concatenate([cf_plate, cf_ramp])
# cf_mask = np.isfinite(cf)

fig, ax = plt.subplots(1,1)
# ax.plot(holden_cp[:,0], 2*holden_cp[:,1], "bo", label="Holden \& Moselle\n(1970, experiment)")
# ax.plot(hung_cp[:,0], 2*hung_cp[:,1], "g^", label="Hung \& MacCormack\n(1976, simulation)")
ax.plot(x_holden, 2*holden_cp[:,1], "bo", label="Holden \& Moselle\n(1970) experiment")
ax.plot(x_hung, 2*hung_cp[:,1], "g^", label="Hung \& MacCormack\n(1976) simulation")
ax.plot(rudy_cp_2d[:,0], 10**(rudy_cp_2d[:,1])/50, "ms", label="Rudy et al (1991)\nsimulation")
ax.plot(x_wall[cp_mask], cp[cp_mask], "r-", label="PLENS")
# ax.set_xlabel(r"$s_{\textrm{wall}}/L$")
ax.set_xlabel(r"$x/L$")
ax.set_ylabel(r"$\displaystyle\frac{2p}{\rho_\infty u_\infty^2}$", rotation=0, labelpad=20)
ax.set_yscale("log")
ax.grid()
ax.legend(loc="best")
fig.tight_layout()
for fmt in ["png", "pdf"]:
    full_figname = "{}cp_comparison_{}.{}".format(res_dir, counter, fmt)
    fig.savefig(full_figname, format=fmt)
    print("Written file {}".format(full_figname))
plt.show()

# fig, ax = plt.subplots(1,1)
# ax.plot(holden_ch[:,0], holden_ch[:,1], "bo", markersize=4, label="Holden \& Moselle\n(1970, experiment)")
# ax.plot(hung_ch[:,0], hung_ch[:,1], "gx", markersize=4, label="Hung \& MacCormack\n(1976, simulation)")
# ax.plot(s_wall[ch_mask], ch[ch_mask], "r-", label="PLENS")
# ax.set_xlabel(r"$s_{\textrm{wall}}/L$")
# ax.set_ylabel(
#     # r"$\displaystyle\frac{-q^{\prime\prime}_y \sec \theta_i}"
#     # r"{\rho_\infty u_\infty \left( h_\infty - h_w + \frac{u_\infty^2}{2} \right)}$",
#     # rotation=0,
#     # labelpad=20
#     r"$c_H$"
# )
# ax.set_yscale("log")
# ax.grid(which="both")
# ax.legend(loc="best")
# fig.tight_layout()
# for fmt in ["png", "pdf"]:
#     full_figname = "{}ch_comparison_{}.{}".format(res_dir, counter, fmt)
#     fig.savefig(full_figname, format=fmt)
#     print("Written file {}".format(full_figname))
# plt.show()

# fig, ax = plt.subplots(1,1)
# ax.plot(holden_cf[:,0], holden_cf[:,1], "bo", markersize=4, label="Holden \& Moselle\n(1970, experiment)")
# ax.plot(hung_cf[:,0], hung_cf[:,1], "gx", markersize=4, label="Hung \& MacCormack\n(1976, simulation)")
# ax.plot(s_wall[cf_mask], cf[cf_mask], "r-", label="PLENS")
# ax.set_xlabel(r"$s_{\textrm{wall}}/L$")
# ax.set_ylabel(
#     # r"$\frac{2\tau_{xy}}{\rho_\infty u_\infty^2}$"
#     # rotation=0,
#     # labelpad=20
#     r"$c_f$"
# )
# ax.grid()
# ax.legend(loc="best")
# fig.tight_layout()
# for fmt in ["png", "pdf"]:
#     full_figname = "{}cf_comparison_{}.{}".format(res_dir, counter, fmt)
#     fig.savefig(full_figname, format=fmt)
#     print("Written file {}".format(full_figname))
# plt.show()
