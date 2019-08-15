#!/usr/bin/env python

import logging
import os
import glob

ROOMT = 298.15
PH2KCAL = 1.364
KCAL2KT = 1.688
KJ2KCAL = 0.239

logging.basicConfig(level=logging.DEBUG, format='%(name)s: %(levelname)s: %(message)s')
logger_env = logging.getLogger("init")
class Env:
    def __init__(self):
        # Hard coded values
        self.runprm = "run.prm"
        self.version = "PyMCCE 1.0"
        self.fn_conflist1 = "head1.lst"
        self.fn_conflist2 = "head2.lst"
        self.fn_conflist3 = "head3.lst"
        self.energy_dir = "energies"

        # load run.prm
        self.prm = self.load_runprm()
        self.prm_default()

        # load ftpl files
        self.tpl = {}
        self.ftpldir = self.prm["TPL_FOLDER"]
        self.load_ftpldir()
        self.tpl_default()
        return

    def load_runprm(self):
        # All values are stripped string
        prm = {}
        logger_env.info("   Loading %s" % self.runprm)
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

        logger_env.info("   Loading ftpl file %s" % file)
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
            if keys in self.tpl:
                if value_string == "!!!":
                # Special rule to negate a key, it is used by the last loaded ftpl file
                # to overwrite values that might have been defined before.
                    del self.tpl[keys]
                else:
                    self.tpl[keys] = self.tpl[keys] + " , " + value_string
            else:
                self.tpl[keys] = value_string
        return


    def prm_default(self):
        if "TPL_FOLDER" not in self.prm or self.prm["TPL_FOLDER"] == "DEFAULT":
            path = str(os.path.dirname(os.path.abspath(__file__)))
            tpl_path = 'param'.join(path.rsplit('bin', 1))
            self.prm["TPL_FOLDER"] = tpl_path
            logger_env.info("   Changing TPL_FOLDER to %s" % tpl_path)

        if "PW_WARNING" not in self.prm:
            self.prm["PW_WARNING"] = "0.05"
            logger_env.info("   Setting PW_WARNING to default value 0.05")

        return


    def tpl_default(self):
        logger_env.info("   Set non-zero default values for missing parameters")
        default_values_keys = [("SCALING", "VDW0"),
                               ("SCALING", "VDW1"),
                               ("SCALING", "VDW"),
                               ("SCALING", "TORS"),
                               ("SCALING", "ELE"),
                               ("SCALING", "DSOLV")]
        for element in default_values_keys:
            if element not in self.tpl:
                logger_env.info("   Set to default: %s = 1.0" % ",".join(element))
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

    def load_ftpldir(self):
        cwd = os.getcwd()
        os.chdir(self.ftpldir)
        files = glob.glob("*.ftpl")
        for fname in files:
            self.load_ftpl(fname)

        os.chdir(cwd)
        return

env = Env()
