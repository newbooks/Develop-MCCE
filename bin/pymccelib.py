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
        self.confnames = {}  # conformer names indexed by residue name
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
        logging.debug("      %s" % file)

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
            if keys[0] == "CONFLIST":
                res = keys[1]
                self.confnames[res] = [x.strip() for x in self.tpl[keys].split(",")]
                for conf in self.confnames[res]:
                    if conf not in self.atomnames:
                        self.atomnames[conf] = []

            elif keys[0] == "CONNECT":
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

    def ftpl2tpl(self, fname):
        """Convert env.tpl to mcce tpl file"""

        # Get a list of residues
        residues = []
        lines = []

        for key in self.tpl:
            if key[0] == "CONFLIST":
                residues.append(key[1])

        for residue in residues:
            # CONFLIST
            conformers = self.confnames[residue]
            value = ["%5s" % x.strip() for x in conformers]
            lines.append("CONFLIST %s        %s\n\n" % (residue, " ".join(value)))

            for conf in self.confnames[residue]:
                # NATOM
                line = "NATOM    %5s      %d\n" % (conf, len(self.atomnames[conf]))
                lines.append(line)
            lines.append("\n")

            # IATOM
            for conf in self.confnames[residue]:
                atoms = self.atomnames[conf]
                iatom = 0
                for atom in atoms:
                    iatom += 1
                    
            # ATOMNAME

            # PROTON

            # PKA

            # ELECTRON

            # EM

            # RXN

            # CONNECT

            # RADIUS

            # CHARGE

            # ROTAMER


        # append 00always
        # scaling factor to tpl record

        open(fname, "w").writelines(lines)
        return

class Atom:
    def __init__(self):
        self.atomname = ""
        self.confname = ""
        self.resname = ""
        self.on = False
        self.iatom = "0"
        return

class Conformer:
    def __init__(self):
        self.confname = ""
        self.resname = ""
        self.atoms = []
        return

class Residue:
    def __init__(self):
        self.resname = []
        self.conformers = []
        return

class Protein:
    """Protein structure"""
    def __init__(self):
        self.residues = []
        return


    def pdb2mcce(self, pdb):
        """Convert pdb to mcce pdb"""
        atom_exceptions = [" H2 ", " OXT", " HXT"]
        mccelines = []
        lines = [x for x in open(pdb).readlines() if x[:6] == "ATOM  " or x[:6] == "HETATM"]

        icount = 0
        previous_resid = ()
        possible_confs = []
        for line in lines:
            # pdb line
            atomname = line[12:16]
            resname = line[17:20]
            chainid = line[21]
            seqnum = int(line[22:26])
            icode = line[26]
            xyz = line[30:54]

            if resname == "PRO" and atomname == " H  ":   # Proline H atom exception
                continue

            current_resid = (resname, chainid, seqnum, icode)
            # mcce line, need to add conf_number, radius, charge, conf_type, conf_history
            if current_resid != previous_resid:
                possible_confs = [x.strip() for x in env.tpl[("CONFLIST", resname)].split(",")]
                logging.info("Identified a new residue %s: %s" % (resname, ", ".join(possible_confs)))
                previous_resid = current_resid
            Found = False
            for confname in possible_confs:
                if atomname in env.atomnames[confname]:
                    conf_type = confname[3:5]
                    conf_number = possible_confs.index(confname)
                    cname = confname
                    Found = True
                    break
            if not Found:
                # this atom is not found in all conformers
                if atomname not in atom_exceptions:
                    print("Atom \"%s\" in pdb file %s can not be assigned to any conformer" % (atomname, pdb))
                continue

            key = ("RADIUS", cname, atomname)
            if key in env.tpl:
                radius_str = env.tpl[key]
                rad, _, _ = radius_str.split(",")
                rad = float(rad)
            else:
                rad = 0.0

            key = ("CHARGE", cname, atomname)
            if key in env.tpl:
                charge_str = env.tpl[key]
                crg = float(charge_str)
            else:
                crg = 0.0

            conf_history = "________"
            newline = "ATOM  %5d %4s %s %c%4d%c%03d%s%8.3f    %8.3f      %s%s\n" % \
                      (icount, atomname, resname, chainid, seqnum, icode, conf_number, xyz, rad, crg, conf_type, conf_history)
            mccelines.append(newline)
            icount += 1

        return mccelines


env = Env()
