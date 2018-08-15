import psutil
import argparse
import datetime
import csv
import time


parser = argparse.ArgumentParser()
parser.add_argument('--file',
                    help='input file name (required)',
                    type=str)
args = parser.parse_args()

# pull args out
filename = args.file

with open(filename) as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    mem_used = 0
    for row in reader:
        mem_row = int(row[4])/1024/1024/1024
        mem_used = max(mem_row,mem_used)

    mem_str =  "Max memory used: %.3f GB" % mem_used
    print(mem_str)
