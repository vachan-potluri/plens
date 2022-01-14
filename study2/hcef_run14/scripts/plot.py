import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams["font.size"] = 14
plt.rcParams["figure.figsize"] = 7,6

p_inf = 47.55513
rho_inf = 0.001206
u_inf = 2667.0
L = 10.17e-2

result_dir = "./group7/trial4/result/"
counter = 9080
cylinder_data = np.genfromtxt(
    "{}/cylinder_data_{}.csv".format(result_dir, counter), delimiter=",", names=True
)
flare_data = np.genfromtxt(
    "{}/flare_data_{}.csv".format(result_dir, counter), delimiter=",", names=True
)
harvey_cp = np.genfromtxt("harvey_run8_cp.csv", delimiter=",")
gnoffo_cp = np.genfromtxt("gnoffo_run8_cp.csv", delimiter=",")
harvey_St = np.genfromtxt("harvey_run8_St.csv", delimiter=",")
gnoffo_St = np.genfromtxt("gnoffo_run8_St.csv", delimiter=",")

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

txy_cylinder = cylinder_data["txy"]
txy_flare = flare_data["txy"]
txy = np.concatenate([txy_cylinder, txy_flare])
txy_mask = np.isfinite(txy)

"""
fig, ax = plt.subplots(1,1)
ax.plot(harvey_cp[:,0], harvey_cp[:,1], "bo", markersize=4, label="Harvey et al (2001)\nexperiment")
ax.plot(gnoffo_cp[:,0], gnoffo_cp[:,1], "gx", markersize=4, label="Gnoffo (2001)\nsimulation")
ax.plot(x[cp_mask], cp[cp_mask], "r-", label="PLENS")
ax.set_xlabel(r"$x/L$")
ax.set_ylabel(r"$\dfrac{p-p_\infty}{\frac{1}{2}\rho_\infty u_\infty^2}$", rotation=0, labelpad=20)
ax.grid()
ax.legend(loc="best")
fig.tight_layout()
fig.savefig("{}/cp_comparison_{}.png".format(result_dir, counter), format="png")
fig.savefig("{}/cp_comparison_{}.pdf".format(result_dir, counter), format="pdf")
plt.show()

del fig, ax
fig, ax = plt.subplots(1,1)
ax.plot(harvey_St[:,0], harvey_St[:,1], "bo", markersize=4, label="Harvey et al (2001)\nexperiment")
ax.plot(gnoffo_St[:,0], gnoffo_St[:,1], "gx", markersize=4, label="Gnoffo (2001)\nsimulation")
ax.plot(x[St_mask], St[St_mask], "r-", label="PLENS")
ax.set_xlabel(r"$x/L$")
ax.set_ylabel(r"$\dfrac{2q^{\prime\prime}_w}{\rho_\infty u_\infty^3}$", rotation=0, labelpad=20)
ax.grid()
ax.legend(loc="best")
fig.tight_layout()
fig.savefig("{}/St_comparison_{}.png".format(result_dir, counter), format="png")
fig.savefig("{}/St_comparison_{}.pdf".format(result_dir, counter), format="pdf")
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
fig.savefig("{}/txy_comparison_{}.png".format(result_dir, counter), format="png")
fig.savefig("{}/txy_comparison_{}.pdf".format(result_dir, counter), format="pdf")
plt.show()
