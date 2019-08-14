#!/usr/bin/env python

import logging

Delta_PW_warning = 0.1
ROOMT = 298.15
PH2KCAL = 1.364
KCAL2KT = 1.688
KJ2KCAL = 0.239

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
        self.tpl = {}
        return

    def load_runprm(self):
        # default is string unless specified as below
        float_values = ["EPSILON_PROT", "TITR_PH0", "TITR_PHD", "TITR_EH0", "TITR_EHD", "CLASH_DISTANCE",
                        "BIG_PAIRWISE", "MONTE_T", "MONTE_REDUCE"]
        int_values = ["TITR_STEPS", "MONTE_RUNS", "MONTE_TRACE", "MONTE_NITER", "MONTE_NEQ",
                      "MONTE_NSTART", "MONTE_FLIPS", "NSTATE_MAX", "MONTE_NEQ"]
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
                    value = fields[0]
                    if key in float_values:
                        prm[key] = float(value)
                    elif key in int_values:
                        prm[key] = int(value)
                    else:
                        prm[key] = value
        return prm

    def print_runprm(self):
        for key in self.prm.keys():
            print("%-25s:%s" % (key, str(self.prm[key])))
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


    def load_tpl(self, file):
        """Load a tpl file."""
        print("   Loading tpl file %s" % file)
        float_values = ["EXTRA", "SCALING"]
        int_values = []

        lines = open(file).readlines()
        for line in lines:
            line = line.split("#")[0]
            if len(line) < 21:
                continue
            keys = [line[:9], line[9:14], line[15:19]]
            value_string = line[20:].strip()

            keys = [x for x in keys if x]
            keys = tuple(keys)

            if keys[0] in float_values:
                self.tpl[keys] = float(value_string)
            elif keys[0] in int_values:
                self.tpl[keys] = int(value_string)
            else:
                self.tpl[keys] = value_string

        return


    def read_extra(self):
        """Read extra.tpl."""
        fname = self.prm["EXTRA"]

        print("   Extra tpl parameters in file %s" % fname)
        if os.path.isfile(fname):
            if fname[-5:] == ".ftpl":
                self.load_ftpl(fname)
            elif fname[-4:] == ".tpl ":
                self.load_tpl(fname)

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