# A python script to delete user-specified result files
# Mainly used for reducing the size of result directory
# Assumes that the file naming is as done by dealii's `write_pvtu_record()` function

import os
import numpy as np

# result directory
res_dir = "/home/vachan/Documents/Work/plens/validation/hypersonic_bl1/N3/trial2/result"

# base name of the output files
output_basename = "output"

# number of processors, and the number of digits used for processor id in output file
n_proc = 64
n_proc_digits = len(str(n_proc))

# extension for archive files
ar_extensions = ["ar", "ar_fixed.data", "ar.info"]

# start and end counters of __all__ the outputs
start_counter = 0
end_counter = 1539

# the counters of outputs to be retained
retain_counters = np.linspace(start_counter, end_counter, 100)
retain_counters = np.append(retain_counters, np.array([end_counter]))

# the counters to be deleted
delete_counters = []

print("Searching for files in {}".format(res_dir))

# search for the files
cur_counter = start_counter
while cur_counter <= end_counter:
    if os.path.isfile("{}/{}_{:06d}.pvtu".format(res_dir, output_basename, cur_counter)):
        if cur_counter in retain_counters:
            print("Found files for counter {}, to be retained".format(cur_counter))
        else:
            print("Found files for counter {}, to be DELETED".format(cur_counter))
            delete_counters.append(cur_counter)
    else:
        print("Not found files for counter {}, ignoring".format(cur_counter))
    cur_counter += 1

# print message and take confirmation
print("The following counters will be deleted")
print(delete_counters)
print("The remaining counters will be retained")

response = input("Proceed with the deletion? [Y/n]: ")
if response == "Y":
    for cur_counter in delete_counters:
        # check again that current counter is not in retain list
        if cur_counter in retain_counters:
            print(
                "Some mistake occurred. " +
                "Counter {} is present in both 'retain' and 'delete' list. ".format(cur_counter) +
                "Skipping this counter."
            )
        else:
            print("Deleting files for counter {}".format(cur_counter))

            # pvtu file
            pvtu_filename = "{}/{}_{:06d}.pvtu".format(res_dir, output_basename, cur_counter)
            if os.path.isfile(pvtu_filename):
                os.remove(pvtu_filename)
            else:
                print("Unable to find file {}".format(pvtu_filename))
            
            # vtu files
            for i in range(n_proc):
                vtu_filename = "{}/{}_{:06d}.{}.vtu".format(
                        res_dir, output_basename, cur_counter, str(i).zfill(n_proc_digits)
                    )
                if os.path.isfile(vtu_filename):
                    os.remove(vtu_filename)
                else:
                    print("Unable to find file {}".format(vtu_filename))
            
            # archive files
            for ext in ar_extensions:
                filename = "{}/{}_{:06d}.{}".format(res_dir, output_basename, cur_counter, ext)
                if os.path.isfile(filename):
                    os.remove(filename)
                else:
                    print("Unable to find file {}".format(filename))
