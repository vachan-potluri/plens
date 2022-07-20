# A script to plot high resolution image of 2d/3d Gmsh mesh (in msh2 format)
# SALOME meshes can be plotted by exporting them as MED (v4.0) files, and opening these files in
# Gmsh to eventually save them in msh2 format (see notable note on SALOME mesh import for FLEXI)

# This was written for high quality vector images of meshes for Manuscript 1 on 20-Jul-2022.

# This script works on a very naive algorithm. Firstly, meshio's mesh object has an attribute called
# "cells" which is a python list containing some or all the following depending on the input msh
# file
# - Physical point indices
# - Physical line connectivity
# - Physical surface connectivity
# - Element connectivity
# The points themselves are stored in the attribute "points". Now we only need the element
# connectivity. In some cases, element connectivity is split into two separate entries in the list.
# To deal with all possible cases, this is the algo followed:
# - Initialise the element connectivity array to the last entry in "cells"
# - Loop over all previous entries
#   - If an entry of the previous list entry's data has the same length as an entry of the element
#     connectivity array (initialised previously)
#     - Add the data of this entry to the element connectivity array

# For 3d meshes, element connectivity has 8 entries we assume the first 4 of them make up the face
# in xy plane.

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

mesh_filename = "/home/vachan/Documents/Work/plens/validation/mengaldo_subsonic_bl/mesh2.msh"
# mesh_filename = "/home/vachan/Documents/Work/plens/validation/jacobs_supersonic_bl2/mesh.msh"
# mesh_filename = "/home/vachan/Documents/Work/plens/study2/lewis_M6_Re1.5e5_adiabatic/lewis_M6_Re1.5e5_adiabatic_trial6.msh"
# mesh_filename = "/home/vachan/Documents/Work/plens/study2/hcef_run9/hcef_run14_trial2.msh"
mesh_file = meshio.read(mesh_filename)
points = mesh_file.points # points array with z coordinate (in general)
points_2d = points[:, :2]

# last item of "cells" will always contain element connectivity
elem_conn = np.array(mesh_file.cells[-1].data)
# now in some cases, some previous items of "cells" can also have element connectivity
# this patch of code is for those cases
if len(mesh_file.cells) > 1:
    # start loop from last but one item of "cells"
    item_id = len(mesh_file.cells)-1
    while item_id > 0:
        # compare the length of an entry in mesh_file.cells[item_id] with an entry from
        # mesh_file.cells[-1]
        cur_cells_item_data = mesh_file.cells[item_id].data
        if len(cur_cells_item_data[0]) == len(elem_conn[0]):
            # this block also has element connectivity
            elem_conn = np.vstack((elem_conn, cur_cells_item_data))
        item_id -= 1

fig, ax = plt.subplots(1, 1, layout="constrained")
elem_id = 0
for quad_ids in elem_conn:
    poly = Polygon(
        np.array([
            points_2d[quad_ids[0]],
            points_2d[quad_ids[1]],
            points_2d[quad_ids[2]],
            points_2d[quad_ids[3]]
        ]),
        edgecolor="black",
        facecolor="gray",
        lw=0.05
    )
    ax.add_patch(poly)
    elem_id += 1
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
