#!/usr/bin/env python

import logging
import os

Delta_PW_warning = 0.1
ROOMT = 298.15
PH2KCAL = 1.364
KCAL2KT = 1.688
KJ2KCAL = 0.239

logging.basicConfig(level=logging.DEBUG, %(levelname)10s: %(message)s')

class Env:
    def __init__(self):
        # Hard coded values
        self.runprm = "run.prm"
        self.version = "PyMCCE 1.0"
        self.fn_conflist1 = "head1.lst"
        self.fn_conflist2 = "head2.lst"
        self.fn_conflist3 = "head3.lst"
        self.energy_dir = "energies"
        self.prm = self.load_runprm()
        if "TPL_FOLDER" not in self.prm or self.prm["TPL_FOLDER"] == "DEFAULT":
            path = str(os.path.dirname(os.path.abspath(__file__)))
            tpl_path = 'param'.join(path.rsplit('bin', 1))
            self.prm["TPL_FOLDER"] = tpl_path
            logging.info("   Changing TPL_FOLDER to %s" % tpl_path)

        self.tpl = {}

        return

    def load_runprm(self):
        # All values are stripped string
        prm = {}
        logging.info("   Loading %s" % self.runprm)
        lines = open(self.runprm).readlines()
        # Sample line: "t        step 1: pre-run, pdb-> mcce pdb                    (DO_PREMCCE)"
        for line in lines:
            line = line.strip()
            line = line.split("#")[0]  # This cuts off everything after #
            left_p = line.rfind("(")
            right_p = line.rfind(")")
            if left_p > 0 and right_p > left_p + 1:
                key = line[left_p + 1:right_p]
                fields = line[:left_p].split()
                if len(fields) >= 1:
                    prm[key] = fields[0]

        return prm

    def print_runprm(self):
        for key in self.prm.keys():
            print("%-25s:%s" % (key, self.prm[key]))
        return

    def load_ftpl(self, file):
        """Load a tpl file."""
        float_values = ["EXTRA", "SCALING"]
        int_values = []

        logging.info("   Loading ftpl file %s" % file)
        lines = open(file).readlines()
        for line in lines:
            line = line.split("#")[0]
            fields = line.split(":")
            if len(fields) != 2:
                continue

            key_string = fields[0].strip()
            keys = key_string.split(",")
            keys = [x.strip().strip("\"") for x in keys]
            keys = [x for x in keys if x]
            keys = tuple(keys)

            value_string = fields[1].strip()
            if keys[0] in float_values:
                self.tpl[keys] = float(value_string)
            elif keys[0] in int_values:
                self.tpl[keys] = int(value_string)
            else:
                self.tpl[keys] = value_string
        return

    def set_default_values(self):
        logging.info("   Set non-zero default values for missing parameters")
        default_values_keys = [("SCALING", "VDW0"),
                               ("SCALING", "VDW1"),
                               ("SCALING", "VDW"),
                               ("SCALING", "TORS"),
                               ("SCALING", "ELE"),
                               ("SCALING", "DSOLV")]
        for element in default_values_keys:
            if element not in self.tpl:
                print("      Set to default: %s = 1.0" % ",".join(element))
                self.tpl[element] = 1.0
        return

    def print_scaling(self):
        """Print scaling factors."""
        # print self.param
        print("      Scaling factors:")
        print("      VDW0  = %.3f" % self.tpl[("SCALING", "VDW0")])
        print("      VDW1  = %.3f" % self.tpl[("SCALING", "VDW1")])
        print("      VDW   = %.3f" % self.tpl[("SCALING", "VDW")])
        print("      TORS  = %.3f" % self.tpl[("SCALING", "TORS")])
        print("      ELE   = %.3f" % self.tpl[("SCALING", "ELE")])
        print("      DSOLV = %.3f" % self.tpl[("SCALING", "DSOLV")])
        return


env = Env()