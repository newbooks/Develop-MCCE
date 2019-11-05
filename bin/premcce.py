#!/usr/bin/env python

import sys
from pymccelib import *
import logging

logging.basicConfig(level=logging.DEBUG, format='%(levelname)-s: %(message)s')

if __name__ == "__main__":
    env.init()
    prot = Protein()
    prot.load_nativepdb(env.prm["INPDB"])
    lines = prot.pdblines()
    sys.stdout.writelines(lines)
