import ibslib_pb as ibslib
import pandas as pd

if __name__ == "__main__":
    twissheader = ibslib.GetTwissHeader("b2_design_lattice_1996.twiss")
    twisstable = ibslib.GetTwissTable("b2_design_lattice_1996.twiss")
    twiss = pd.DataFrame.from_dict(twisstable)
    print(ibslib.clight)
    print(ibslib.proton_radius)
