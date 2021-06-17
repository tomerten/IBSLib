import ibslib_pb as ibslib

if __name__ == "__main__":
    twissheader = ibslib.GetTwissHeader("b2_design_lattice_1996.twiss")
    print(twissheader)
