#!/usr/bin/env python

import logging
import os
import glob

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
        self.ftpldir = ""
        self.prm = {}
        self.tpl = {}
        self.atomnames = {}  # atom names indexed by conformer name
        return

    def init(self):
        logging.info("Step 0. Initialize MCCE run environment.")

        # load run.prm
        self.prm = self.load_runprm()
        self.prm_default()

        # load ftpl files
        self.ftpldir = self.prm["TPL_FOLDER"]
        self.load_ftpldir()

        # load extra.ftpl
        if "EXTRA" in self.prm and os.path.isfile(self.prm["EXTRA"]):
            self.load_ftpl(self.prm["EXTRA"])

        # revise self.atomnames to include empty conformer types
        for res_conflist in [x for x in self.tpl.keys() if x[0] == "CONFLIST"]:
            for conf in [x.strip() for x in self.tpl[res_conflist].strip().split(",")]:
                if conf not in self.atomnames:
                    self.atomnames[conf] = []

        logging.info("Step 0. Done.\n")
        return

    def load_runprm(self):
        # All values are stripped string
        prm = {}
        path = os.path.abspath(self.runprm)
        logging.info("   Loading run parameter from %s" % path)
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
        logging.debug("   Loading from file %s" % file)

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
            if keys in self.tpl:   # Overwrite
                logging.warning("   Value of \"%s\": (%s) is replaced by (%s)" % (",".join(keys), self.tpl[keys],
                                                                                   value_string))

            self.tpl[keys] = value_string

            # Make an atom list in the natural order of CONNECT record.
            if keys[0] == "CONNECT":
                atom = keys[1]
                conf = keys[2]
                if conf in self.atomnames:
                    self.atomnames[conf].append(atom)
                else:
                    self.atomnames[conf] = [atom]


        return

    def prm_default(self):
        if "TPL_FOLDER" not in self.prm or self.prm["TPL_FOLDER"].upper() == "DEFAULT":
            path = str(os.path.dirname(os.path.abspath(__file__)))
            tpl_path = 'param'.join(path.rsplit('bin', 1))
            self.prm["TPL_FOLDER"] = tpl_path
            logging.info("   Default TPL_FOLDER is set to %s" % tpl_path)
        if "DELPHI_EXE" not in self.prm or self.prm["DELPHI_EXE"].upper() == "DEFAULT":
            path = str(os.path.dirname(os.path.abspath(__file__)))
            self.prm["DELPHI_EXE"] = path
            logging.info("   Default DELPHI_EXE is set to %s" % path)
        if "SCALING_VDW0" not in self.prm or self.prm["SCALING_VDW0"].upper() == "DEFAULT":
            logging.info("      Set to default: SCALING_VDW0 = 1.0")
            self.prm["SCALING_VDW0"] = "1.0"
        if "SCALING_VDW1" not in self.prm or self.prm["SCALING_VDW1"].upper() == "DEFAULT":
            logging.info("      Set to default: SCALING_VDW1 = 1.0")
            self.prm["SCALING_VDW1"] = "1.0"
        if "SCALING_VDW" not in self.prm or self.prm["SCALING_VDW"].upper() == "DEFAULT":
            logging.info("      Set to default: SCALING_VDW = 1.0")
            self.prm["SCALING_VDW"] = "1.0"
        if "SCALING_TORS" not in self.prm or self.prm["SCALING_TORS"].upper() == "DEFAULT":
            logging.info("      Set to default: SCALING_TORS = 1.0")
            self.prm["SCALING_TORS"] = "1.0"
        if "SCALING_ELE" not in self.prm or self.prm["SCALING_ELE"].upper() == "DEFAULT":
            logging.info("      Set to default: SCALING_ELE = 1.0")
            self.prm["SCALING_ELE"] = "1.0"
        if "SCALING_DSOLV" not in self.prm or self.prm["SCALING_DSOLV"].upper() == "DEFAULT":
            logging.info("      Set to default: SCALING_DSOLV = 1.0")
            self.prm["SCALING_DSOLV"] = "1.0"

        return

    def print_scaling(self):
        """Print scaling factors."""
        # print self.param
        print("      Scaling factors:")
        print("      VDW0  = %.3f" % self.prm["SCALING_VDW0"])
        print("      VDW1  = %.3f" % self.prm["SCALING_VDW1"])
        print("      VDW   = %.3f" % self.prm["SCALING_VDW"])
        print("      TORS  = %.3f" % self.prm["SCALING_TORS"])
        print("      ELE   = %.3f" % self.prm["SCALING_ELE"])
        print("      DSOLV = %.3f" % self.prm["SCALING_DSOLV"])
        return

    def load_ftpldir(self):
        cwd = os.getcwd()
        os.chdir(self.ftpldir)
        files = glob.glob("*.ftpl")
        files.sort()
        logging.info("   Loading ftpl files from %s" % self.ftpldir)
        for fname in files:
            self.load_ftpl(fname)

        os.chdir(cwd)
        return

class Atom:
    def __init__(self):
        self.name = ""
        self.confname = ""
        self.resname = ""
        self.on = False
        self.iatom = "0"
        self.chainid = ""
        self.seqnum = 0
        self.icode = ""
        self.xyz=()
        self.resid = ()
        return

    def load_nativeline(self, line):
        self.name = line[12:16]
        self.resname = line[17:20]
        self.chainid = line[21]
        self.seqnum = int(line[22:26])
        self.icode = line[26]
        self.xyz = (float(line[30:38]), float(line[38:46]),float(line[46:54]))

        self.resid = (self.resname, self.chainid, self.seqnum, self.icode)
        return

class Conformer:
    def __init__(self):
        self.confname = ""
        self.resname = ""
        self.atoms = []
        self.resid = ()
        return

class Residue:
    def __init__(self, resid):
        self.resid = resid
        self.resname = resid[0]
        self.conformers = []
        return

class Protein:
    """Protein structure"""
    def __init__(self):
        self.residues = []
        return


    def load_nativepdb(self, pdb):
        """Load native pdb file."""

        lines = [x for x in open(pdb).readlines() if x[:6] == "ATOM  " or x[:6] == "HETATM"]

        resids = []
        for line in lines:
            # pdb line
            atom = Atom()
            atom.load_nativeline(line)

            # mcce line, need to add conf_number, radius, charge, conf_type, conf_history
            try:
                ires = resids.index(atom.resid)
            except ValueError:  # new residue
                self.residues.append(Residue(atom.resid))
                resids.append(atom.resid)
                ires = len(self.residues) - 1
                self.residues[ires].conformers=[Conformer()]   # set up conformer 0
            # load atoms into conformer 0
            self.residues[ires].conformers[0].atoms.append(atom)

        # separate side chain atoms from backbone - BK atoms remain in conformer 0, the rest go to conformer 1
        for res in self.residues:
            conflist = [x.strip() for x in env.tpl[("CONFLIST", res.resname)].strip().split(",")]
            if res.conformers:
                for atom in res.conformers[0].atoms:
                    # find the first conformer type this atom fits
                    for conftype in conflist:
                        if atom.name in env.atomnames[conftype]:
                            if conftype[-2:] == "BK":  # stays in this conformer, break search conf, next atom
                                break
                            else:
                                if len(res.conformers) > 1:
                                    res.conformers[1].atoms.append(atom)
                                else:
                                    res.conformers.append(Conformer())
                                    res.conformers[1].confname = conftype
                                    res.conformers[1].resname = res.resname
                                    res.conformers[1].atoms.append(atom)

        # delete atoms don't belong to conformer 1
        for res in self.residues:
            if len(res.conformers) > 1:
                confname = res.conformers[1].confname
                valid_atoms = env.atomnames[confname]
                for atom in res.conformers[1].atoms:
                    if atom.name not in valid_atoms:
                        logging.WARNING("   Deleted atom \"%s\" of %s because it doesn't fit into initial conformer." % (
                            atom.name, res.resname))

        return

    def pdblines(self):
        lines = []
        icount = 1
        for res in self.residues:
            conflist = env.tpl[("CONFLIST", res.resname)]
            logging.debug(conflist)
            for conf in res.conformers:
                for atom in conf.atoms:
                    line = "ATOM  %5d %4s %3s %c\n" % (icount, atom.name, atom.resname, atom.chainid)
                    lines.append(line)
                    icount += 1
        return lines


env = Env()
