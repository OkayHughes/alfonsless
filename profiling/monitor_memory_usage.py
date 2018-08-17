pimport psutil
import argparse
import datetime
import csv
import time


def monitor_memory_usage(pid, sleep_time=None, filename = None):

    # pull args out
    p = psutil.Process(pid)

    if sleep_time is None:
        sleep_time = 60 ;

    if filename is None:
        now = datetime.datetime.now()
        year = str(now.year)
        month = '{:02d}'.format(now.month)
        day = '{:02d}'.format(now.day)
        filename = str(pid) + '_' + year + month + day + '.csv'

    t = 1
    while t:
        mem = p.memory_info()
        with open(filename, 'a', newline='') as csvfile:
            now = datetime.datetime.now()
            writer = csv.writer(csvfile)
            writer.writerow([pid, now.hour, now.minute, now.second, mem.rss,
                             mem.vms])
        time.sleep(sleep_time)
        if p.status() == 'zombie':
            t = 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--pid',
                        help='PID of the process that is monitored. (required)',
                        type=int)
    parser.add_argument('-t',
                        help='Sleep time between mem checks in [s], default 60s',
                        type=int,
                        default=60)
    parser.add_argument('-o',
                        help='Output file name',
                        type=str)
    args = parser.parse_args()
    sleep_time = args.t
    filename = args.o
    pid = args.pid

    monitor_memory_usage(pid, sleep_time, filename)

