#!/usr/bin/ipython
from profiling.coll_output import out_call


def graph_node(num_poly, deg, num_consts, pref):
    out_call([num_poly, deg, num_consts], "test_framework", "{}_dv={}_deg={}_nc={}.log".format(pref, num_poly, deg, num_consts))

if __name__ == "__main__":
    from sys import argv
    graph_node(argv[1], argv[2], argv[3], argv[4])
