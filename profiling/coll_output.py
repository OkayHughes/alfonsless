from rec_mat import call_matlab
def handle(fname, val):
    with open(fname, "a") as f:
        f.write(val)

def out_call(params, scname, fname):
    hand = lambda x: handle(fname, x)
    with open(fname, "w") as f:
        pass
    st_params = ", ".join(map(str, params))
    call_matlab("{}({})".format(scname, st_params), hand)

def save_out_call(deg, fname):
    hand = lambda x: handle(fname, x)
    with open(fname, "w") as f:
        pass
    call_matlab("test_out({})".format(deg), hand)

def save_out_call_defective(deg, fname):
        hand = lambda x: handle(fname, x)
        with open(fname, "w") as f:
            pass
        call_matlab("test_out_defective({})".format(deg), hand)

def save_out_call_more_def(deg, fname):
        hand = lambda x: handle(fname, x)
        with open(fname, "w") as f:
            pass
        call_matlab("test_out_more_def({})".format(deg), hand)

