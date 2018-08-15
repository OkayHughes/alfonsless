from profiling.monitor_memory_usage import monitor_memory_usage
from subprocess import Popen
from subprocess import PIPE
from subprocess import STDOUT
import sys
import argparse


def comp_degree(d):
    print("Beginning okay")
    command = "/Applications/MATLAB_R2018a.app/bin/matlab -nosplash -nodisplay -nodesktop -r \"addpath('segway'); segway_FRS_solver_okay({}); exit;\""

    proc = Popen(command.format(d), shell=True, stdin=PIPE, stdout=sys.stdout, stderr=sys.stdout, close_fds=True)

    monitor_memory_usage(proc.pid, 1, "mem_okay_{}.log".format(d));

    print("End okay")
    print("Beginning shreyas")

    command = "/Applications/MATLAB_R2018a.app/bin/matlab -nosplash -nodisplay -nodesktop -r \"addpath('segway'); segway_FRS_solver_shrey({}); exit;\""

    proc = Popen(command.format(d), shell=True, stdin=PIPE, stdout=sys.stdout, stderr=sys.stdout, close_fds=True)
    
    monitor_memory_usage(proc.pid, 1, "mem_shrey_{}.log".format(d));

    print("End shreyas")



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('d',
                        help='Max degree of solution',
                        type=int)
    args = parser.parse_args()

    comp_degree(args.d)