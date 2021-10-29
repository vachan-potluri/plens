# A script to generate a pvd file from existing pvtu files in the current directory.
# Handy when only specific outputs are copied from a remote machine (on which simulation is
# running) for viewing results on a local machine.
# Assumes that the pvtu files are named as "<basename>_<6 digit counter>.pvtu"

import os
import numpy as np

# tab string
tab = "    "

# Returns the formatted string for an entry
def get_entry_string(time, filename):
    temp = '\t\t<DataSet timestep="{}" group="" part="0" file="{}"/>\n'.format(time, filename)
    return temp.expandtabs(4)

pvd_header = (
    '<?xml version="1.0"?>\n'
    '<VTKFile type="Collection" version="0.1" ByteOrder="LittleEndian">\n'
    '\t<Collection>\n'
).expandtabs(4)

pvd_footer = (
    '\t</Collection>\n'
    '</VTKFile>\n'
).expandtabs(4)

counter_list = np.array([], dtype=int)
filename_list = np.array([], dtype=str)

for filename in os.listdir("."):
    if filename.endswith(".pvtu"):
        # get the counter as a string
        counter_str = filename[-11:-5]
        counter = int(counter_str)
        counter_list = np.append(counter_list, [counter])
        filename_list = np.append(filename_list, [filename])

indices = counter_list.argsort()
print("Found counters\n{}".format(counter_list))

pvd_file = open("dummy.pvd", "w")
pvd_file.write(pvd_header)
for i in indices:
    pvd_file.write(get_entry_string(counter_list[i], filename_list[i]))
pvd_file.write(pvd_footer)
pvd_file.close()

print("Written file dummy.pvd")
