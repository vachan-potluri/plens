# A file to delete a fraction of su2 result files
# This script has to be called from the case directory where the config file etc. are accessible
# This script doesn't delete the directories, but generates a bash file which can be executed later
# to delete the directories

import glob
import numpy as np
import os
import subprocess

bash_file = open("delete_su2_files.sh", "w")
file_groups = ["restart", "surface", "vol_solution"]
file_extensions = ["dat", "vtu", "vtu"]

for group,ext in zip(file_groups,file_extensions):
    print(f"Looking for {group} files")
    file_list = glob.glob(f"{group}*{ext}")

    if len(file_list) > 0:
        times = np.array([], dtype=int)
        for file in file_list:
            times = np.append(times, int((file.split(".")[0]).split("_")[-1]))
        times = np.sort(times)
        if times[0] == 0: times = times[1:]
        print("Found time counters:")
        print((", ".join(["{:05d}"]*len(times))).format(*times))
        n_times = len(times)
        retain_indices = np.arange(n_times-1, 0, -10)
        del_indices = list(set(range(n_times)) - set(retain_indices))
        print("Deleting the following time counters")
        print((", ".join(["{:05d}"]*len(del_indices))).format(*times[del_indices]))
        for t_del in times[del_indices]:
            bash_file.write("rm {}_{:05d}.{}\n".format(group,t_del,ext))

print(f"Written delete commands into {bash_file.name}")