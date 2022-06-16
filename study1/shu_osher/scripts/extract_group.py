# Extracts line data for all cases
# Introduced on 16-Jun-2022 for SIS2022 analysis

# This must be executed from the "run" directory from where individual case directories are
# accessible.

import subprocess

case_suffix = "chandrashekhar"
res_dir = "result_16Jun2022"
script_dir = "/home/vachan/Documents/Work/plens/study1/shu_osher/scripts"

for N in [1,2,3,5]:
    for dof in [200,400,800]:
        res_dir_relative = "dof{0}_{1}_{1}_N{2}_{3}/{4}".format(
            dof, N+1, N, case_suffix, res_dir
        ) # current result directory relative to "run" directory
        times_file = res_dir_relative + "/output.times"
        # get the last line
        runlog = subprocess.run(
            ["tail", "-n", "1", times_file],
            stdout=subprocess.PIPE
        ).stdout.decode("utf-8")
        line_splits = runlog.split()
        end_ctr = int(line_splits[0])
        end_time = float(line_splits[1])
        print("In {}, found end counter {}".format(res_dir_relative, end_ctr))
        if end_time != 1.8:
            # prints in yellow
            print(
                "\033[93m" +
                "Warning: the found end time value {}".format(end_time) +
                " is not equal to the required end time 1.8" +
                "\033[0m"
            )
        subprocess.run([
            "pvpython",
            script_dir + "/extract_data.py",
            res_dir_relative,
            "{}".format(end_ctr)
        ])