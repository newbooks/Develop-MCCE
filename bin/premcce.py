#!/usr/bin/env python

import sys
from pymccelib import *
import logging

logging.basicConfig(level=logging.DEBUG, format='%(levelname)-s: %(message)s')

if __name__ == "__main__":
    env.init()
    prot = Protein()
    prot.load_nativepdb(env.prm["INPDB"])

    # identify N and C terminal
    if env.prm["TERMINALS"].upper() == "T":
        prot.identify_nc()

    # remove exposed water

    # Disulfide bridge



    lines = prot.pdblines()
    open(env.fn_step1_out,"w").writelines(lines)
