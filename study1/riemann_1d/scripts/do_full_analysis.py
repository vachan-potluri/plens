# A single script to do the following tasks in the current directory (which is assumed to have all
# the result files):
# - Note the end time and counter from `output.times`
#   - In case simulation had crashed midway, the end time would let us know
# - Note the CPU time and time steps from `plens.log`
# - Extract line data (outsourced)
# - Compute 3d velocity error
# - Plot line data, save the plotted data for group plotting
# - Compute the errors wrt exact solution

# All this data is saved in `full_analysis.log`



import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import re
import subprocess
from ic_data import end_times as test_end_times
from ic_data import test_options



parser = argparse.ArgumentParser()
parser.add_argument(
    "test",
    help="The test case analysed. Options are: {}.".format(test_options)
)
args = parser.parse_args()



def check_file_existance(filename):
    # asserts if a `filename` exists
    assert os.path.isfile(filename), "Unable to find file '{}'".format(filename)

def parse_value_in_string(s, substr, size=10):
    # searches for `substr` in `s` and returns a string of the content that follows this and
    # halts at "," or "\n"
    # eg: given s = "xyz b: 35 , c: ab", calling this with "b" as `substr` returns "35"
    # the length of the value following colon is given by the `size` parameter
    lookstr = substr
    loc = s.find(lookstr)
    assert loc != -1, "Not found substring"
    temp = s[loc+len(lookstr) : loc+len(lookstr)+size]
    temp_splits = re.split(r"[,\n]+", temp)
    return temp_splits[0]



logfile = open("full_analysis.log", "w")
print("Doing full analysis in current directory")
# directory where outsourced scripts lie
script_dir = "/home/vachan/Documents/Work/plens/study1/riemann_1d/scripts/"



## 1. Get the end time and counter from `output.times`
# check if the times file exists
print("\nStep 1")
check_file_existance("output.times")
# get the last line
runlog = subprocess.run(
    ["tail", "-n", "1", "output.times"],
    stdout=subprocess.PIPE
).stdout.decode("utf-8")
line_splits = runlog.split()
end_ctr = int(line_splits[0])
end_time = float(line_splits[1])
print("Found end time {} and end counter {}".format(end_time, end_ctr))
logfile.write("end time, {}\nend counter, {}\n".format(end_time, end_ctr))
if end_time != test_end_times[args.test]:
    # prints in yellow
    print(
        "\033[93m" +
        "Warning: the found end time value {} is not equal to the required end time {}".format(
            end_time,
            test_end_times[args.test]
        ) +
        "\033[0m"
    )



## 2. Get time steps and CPU time from `plens.log`
# check if the log file exists
print("\nStep 2")
check_file_existance("plens.log")
# get the last 30 lines (will contain only the last time step)
runlog = subprocess.run(
    ["tail", "-n", "30", "plens.log"],
    stdout=subprocess.PIPE
).stdout.decode("utf-8")
n_timesteps = int(parse_value_in_string(runlog, "time steps: "))
cpu_time = float(parse_value_in_string(runlog, "CPU time: "))
print("Found # time steps {} and CPU time {}".format(n_timesteps, cpu_time))
logfile.write("time steps, {}\ncpu time, {}\n".format(n_timesteps, cpu_time))
logfile.write("cpu time per time step, {}\n".format(cpu_time/n_timesteps))



## 3. Extract line data using `extract_data.py`
print("\nStep 3")
print("Extracting data (outsourced)")
subprocess.run([
    "pvpython",
    script_dir + "extract_data.py",
    ".",
    str(end_ctr),
    "-r",
    "6400"
])



## 4. Compute 3d velocity error
print("\nStep 4")
print("Computing 3d velocity error (outsourced)")
runlog = subprocess.run(
    ["pvpython", script_dir + "compute_3d_error.py", ".", str(end_ctr)],
    stdout=subprocess.PIPE
).stdout.decode("utf-8")
print(runlog)
abs_vel_error_3d = float(parse_value_in_string(runlog, "error: ", 30))
logfile.write("3d absolute velocity error, {}\n".format(abs_vel_error_3d))



## 5. Plotting line data and comparing with exact solution
# first getting diaphragm location from IC (0th output)
print("\nStep 4")
jump_vars = {
    "test1-1": "p",
    "test1-2": "u",
    "test1-3": "p",
    "test1-4": "p",
    "test1-5": "p"
}
# get line data for ic
print("Getting line data for IC to detect diaphragm position (outsourced)")
subprocess.run(
    ["pvpython", script_dir + "extract_data.py", ".", "0", "-r", "6400"]
)
ic_data = np.genfromtxt("line_data_000000.csv", delimiter=",", names=True)
x = ic_data["Points0"]
ic_vec = ic_data[jump_vars[args.test]]
jump_tol = 1e-4
jump_indices = []
i=0
for i in range(len(ic_vec)-1):
    if abs(ic_vec[i+1]-ic_vec[i]) > abs(ic_vec[i])*jump_tol:
        jump_indices.append(i)
print("Found jumps in variable {} at locations:".format(jump_vars[args.test]))
print(x[jump_indices])
plt.plot(x, ic_vec, "bo", markersize=1)
plt.xlabel(r"$x$")
plt.ylabel(jump_vars[args.test])
plt.title("IC")
plt.grid()
plt.show()
dia_loc = x[jump_indices][0]
print("Using diaphragm location {}".format(dia_loc))
runlog = subprocess.run(
    [
        "python3",
        script_dir + "plot_individual.py",
        "./line_data_{:06d}.csv".format(end_ctr), # "./" is required to detect dir for saving fig
        args.test,
        str(dia_loc),
        "comparison"
    ],
    stdout=subprocess.PIPE
).stdout.decode("utf-8")
print(runlog)
l1_error = float(parse_value_in_string(runlog, "1, "))
l2_error = float(parse_value_in_string(runlog, "2, "))
linf_error = float(parse_value_in_string(runlog, "inf, "))
logfile.write(
    "l1 error, {}\nl2 error, {}\nlinf error, {}\n".format(l1_error, l2_error, linf_error)
)



logfile.close()
