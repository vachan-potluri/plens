# A file to delete a fraction of openfoam result files
# This script has to be called from the case directory where folders 0, system and constant are
# accessible
# This script doesn't delete the directories, but generates a bash file which can be executed later
# to delete the directories

# Keep in mind that openfoam log file takes the largest space

import glob
import numpy as np
import os
import subprocess

def delete_time_dirs(file, times, n_del):
    # this function deletes selective time directories from the list "times" in the current directory
    # it is checked that "times" doesn't contain "0"
    # only one in every "n_del" time directories is retained
    # the last time directory is not deleted
    if times[0] == "0":
        times = times[1:]
    n_times = len(times)

    # return if only one non-zero time directory is seen
    if n_times == 1: return

    retain_indices = np.arange(n_times-1,0,-n_del)
    del_indices = list(set(range(n_times)) - set(retain_indices))
    print("Deleting the following time directories:")
    print(times[del_indices])
    for t_del in times[del_indices]:
        file.write(f"rm -r {t_del}\n")


# first delete the reconstructed (serialised) time directories
bash_file = open("delete_openfoam_files.sh", "w")
subprocess_out = subprocess.run(
    ["foamListTimes"],
    stdout = subprocess.PIPE,
    text=True
)
serial_times = np.array(subprocess_out.stdout.split("\n")[:-1])
print("Found {} serial times:".format(len(serial_times)), serial_times)
delete_time_dirs(bash_file, serial_times, 10)

# check if this is a parallel run
processor_dirs = glob.glob("processor*/")
if len(processor_dirs) != 0:
    print("Detected a parallel run")
    print("Finding number of processors ...")
    n_proc = 0
    for proc_dir in processor_dirs:
        n_proc = max(int((proc_dir.split("/")[0]).split("processor")[1]), n_proc)
    print("\tParallel run performed with {} processors".format(n_proc+1))
    subprocess_out = subprocess.run(
        ["foamListTimes", "-processor"],
        stdout = subprocess.PIPE,
        text=True
    )
    parallel_times = np.array(subprocess_out.stdout.split("\n")[:-1])
    print("Found {} parallel times:".format(len(parallel_times)), parallel_times)

    for i in range(n_proc+1):
        cmd = "cd processor{}".format(i)
        bash_file.write(f"{cmd}\n")
        subprocess.run(cmd, shell=True)
        print("In directory processor{}".format(i))
        delete_time_dirs(bash_file, parallel_times, 10)
        cmd = "cd .."
        bash_file.write(f"{cmd}\n")
        subprocess.run(cmd, shell=True)

# look for log file
log_files = glob.glob("*log*")
if len(log_files) != 0:
    print(f"Found the following log files:\n{log_files}")
    del_log_file = input("Do you want to delete the log files? [y/n]: ")
    if del_log_file == "y":
        bash_file.write("rm {}\n".format(" ".join(log_files)))
print(f"Written delete commands into {bash_file.name}")