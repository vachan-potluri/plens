# A script to plot high resolution 2d image of 2d/3d Gmsh mesh (in msh2 format). This is designed
# for plotting 2d meshes with or without extrusion in the z-direction.
# SALOME meshes can be plotted by exporting them as MED (v4.0) files, and opening these files in
# Gmsh to eventually save them in msh2 format (see notable note on SALOME mesh import for FLEXI)

# This was written for high quality vector images of meshes for Manuscript 1 on 20-Jul-2022.

# This script loops over all items of the "cells" attribute which are of type "quad", and draws
# them. This may cause some issues for axisymmetric meshes that are not pure z-projections.

import meshio
import matplotlib.pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
    "font.size": 10,
    "axes.formatter.limits": [-2,2]
})
from matplotlib.patches import Polygon
import numpy as np

# mesh_filename = "/home/vachan/Documents/Work/plens/validation/mengaldo_subsonic_bl/mesh2.msh"
# mesh_filename = "/home/vachan/Documents/Work/plens/validation/jacobs_supersonic_bl2/mesh.msh"
mesh_filename = "/home/vachan/Documents/Work/plens/study2/lewis_M6_Re1.5e5_adiabatic/lewis_M6_Re1.5e5_adiabatic_trial6.msh"
# mesh_filename = "/home/vachan/Documents/Work/plens/study2/hcef_run9/hcef_run14_trial2.msh"
mesh_file = meshio.read(mesh_filename)
points = mesh_file.points # points array with z coordinate (in general)
points_2d = points[:, :2]

# initialise the quad connectivity array to be used for plotting
quad_conn = np.zeros((1,4), dtype=int)
# loop over all items in the "cells" list and look for an item of type "quad"
for c in mesh_file.cells:
    if c.type == "quad":
        quad_conn = np.vstack((quad_conn, c.data))
# delete the first row that was set fot initialisation
quad_conn = np.delete(quad_conn, (0), axis=0)

fig, ax = plt.subplots(1, 1, layout="constrained")
for qid, quad_vertex_ids in enumerate(quad_conn):
    # if np.all(np.isclose(points[quad_vertex_ids][2], points[quad_vertex_ids[0]][2])) == False:
    #     print("WARNING: quad {} has some points out of xy plane".format(qid))
    poly = Polygon(
        np.array([
            points_2d[quad_vertex_ids[0]],
            points_2d[quad_vertex_ids[1]],
            points_2d[quad_vertex_ids[2]],
            points_2d[quad_vertex_ids[3]]
        ]),
        edgecolor="black",
        facecolor="gray",
        lw=0.05
    )
    ax.add_patch(poly)
ax.set_aspect("equal")
ax.set_xlim(np.min(points_2d[:,0]), np.max(points_2d[:,0]))
ax.set_ylim(np.min(points_2d[:,1]), np.max(points_2d[:,1]))
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$", rotation=0, labelpad=5)
# https://stackoverflow.com/a/73047823/8028836
fig.draw_without_rendering()
tb = fig.get_tightbbox(fig.canvas.get_renderer())
fig.set_size_inches(tb.width, tb.height)
# fig.tight_layout(pad=0.25)
fig.savefig("mesh.pdf", format="pdf")
print("(Over)Written file mesh.pdf")
plt.show()
