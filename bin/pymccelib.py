#!/usr/bin/env python

import logging
import os
import glob
from geometry import *

ROOMT = 298.15
PH2KCAL = 1.364
KCAL2KT = 1.688
KJ2KCAL = 0.239
DEFAULT_RAD = 1.5  # for dielectric boundary
AMINO_ACIDS = ["ALA", "ARG", "ASN", "ASP", "CYS", "CYL", "GLN", "GLY", "GLU", "HIS", "HIL", "ILE", "LEU", "LYS",
               "MET", "PHE", "PRO", "THR", "TRP", "TYR", "VAL"]



class Env:
    def __init__(self):
        # Hard coded values
        self.runprm = "run.prm"
        self.version = "PyMCCE 1.0"
        self.fn_conflist1 = "head1.lst"
        self.fn_step1_out = "step1_out.pdb"
        self.fn_conflist2 = "head2.lst"
        self.fn_step2_out = "step2_out.pdb"
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
        # Sample line: "t   include     step 1: pre-run, pdb-> mcce pdb                    (DO_PREMCCE)"
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
        self.icount = 0
        self.iconf = 0
        self.name = ""
        self.confname = ""
        self.resname = ""
        self.on = False
        self.iatom = "0"
        self.chainid = ""
        self.seqnum = 0
        self.icode = "_"
        self.xyz=()
        self.resid = ()
        self.crg = 0.0
        self.rad = DEFAULT_RAD
        self.history = "__________"
        return

    def load_nativeline(self, line):
        self.name = line[12:16]
        self.resname = line[17:20]
        self.chainid = line[21]
        self.seqnum = int(line[22:26])
        if line[26] != " ":
            self.icode = line[26]
        self.xyz = (float(line[30:38]), float(line[38:46]),float(line[46:54]))

        self.resid = (self.resname, self.chainid, self.seqnum, self.icode)
        return

    def printme(self):
        line = "ATOM  %5d %4s %3s %c%04d%c%03d%8.3f%8.3f%8.3f%8.3f    %8.3f      %s\n" % (self.icount, self.name,
                                                                             self.resname,
                                                                     self.chainid,
                                                 self.seqnum, self.icode, self.iconf, self.xyz[0],
                                                 self.xyz[1], self.xyz[2], self.rad, self.crg, self.history)

        return line

class Conformer:
    def __init__(self):
        self.confname = ""
        self.resname = ""
        self.atoms = []
        self.resid = ()
        self.history = "__________"
        return

class Residue:
    def __init__(self, resid):
        self.resid = resid
        self.resname, self.chainid, self.seqnum, self.icode = resid
        conf = Conformer()
        conf.history = "BK________"
        self.conformers = [conf]
        self.flag = ""    # flag for ntr, ctr label or other purpose
        return

class Protein:
    """Protein structure"""
    def __init__(self):
        self.residues = []
        return


    def load_nativepdb(self, pdb):
        """Load native pdb file."""

        lines = [x for x in open(pdb).readlines() if x[:6] == "ATOM  " or x[:6] == "HETATM"]

        atoms = []
        for line in lines:
            # pdb line
            atom = Atom()
            atom.load_nativeline(line)
            atoms.append(atom)

        self.load_atoms_single(atoms)

        return

    def load_atoms_single(self, atoms):
        resids = []
        self.residues = []
        for atom in atoms:
            try:
                ires = resids.index(atom.resid)
            except ValueError:  # new residue
                self.residues.append(Residue(atom.resid))
                resids.append(atom.resid)
                ires = len(self.residues) - 1

            # load atoms into conformer 0
            self.residues[ires].conformers[0].atoms.append(atom)

        # separate side chain atoms from backbone - BK atoms remain in conformer 0, the rest go to conformer 1
        for res in self.residues:
            conflist = [x.strip() for x in env.tpl[("CONFLIST", res.resname)].strip().split(",")]
            if res.conformers:
                new_conf0 = []
                for atom in res.conformers[0].atoms:
                    # find the first conformer type this atom fits
                    for conftype in conflist:
                        if atom.name in env.atomnames[conftype]:
                            if conftype[-2:] == "BK":  # stays in this conformer, break search conf, next atom
                                new_conf0.append(atom)
                            else:
                                if len(res.conformers) > 1:
                                    res.conformers[1].atoms.append(atom)
                                else:
                                    conf = Conformer()
                                    conf.history = "%2s________" % (conftype[-2:])  # last two characters
                                    res.conformers.append(conf)
                                    res.conformers[1].confname = conftype
                                    res.conformers[1].resname = res.resname
                                    res.conformers[1].atoms.append(atom)
                            break  # do not search other conformers
                res.conformers[0].atoms = new_conf0


        # delete atoms don't belong to conformer 1
        for res in self.residues:
            if len(res.conformers) > 1:
                confname = res.conformers[1].confname
                valid_atoms = env.atomnames[confname]
                conf1_atoms = []
                for atom in res.conformers[1].atoms:
                    if atom.name in valid_atoms:
                        conf1_atoms.append(atom)
                    else:
                        logging.WARNING("   Deleted atom \"%s\" of %s because it doesn't fit into initial conformer." % (
                            atom.name, res.resname))

        return

    def pdblines(self):
        lines = []
        icount = 1
        for res in self.residues:
            conflist = env.tpl[("CONFLIST", res.resname)]
            #logging.debug(conflist)
            iconf = 0
            for conf in res.conformers:
                for atom in conf.atoms:
                    atom.icount = icount
                    atom.iconf = iconf
                    atom.history = conf.history
                    line = atom.printme()
                    lines.append(line)
                    icount += 1
                iconf += 1
        return lines

    def identify_nc(self):
        # The first and last amino acid in a chain, and no extra bonded atom from other residues in the same chain
        clash_distance = float(env.prm["CLASH_DISTANCE"])
        clash_distance2 = clash_distance * clash_distance

        confnames = [x.strip() for x in env.tpl[("CONFLIST", "NTR")].split(",")]
        NTR_atomnames = set()
        for conf in confnames:
            NTR_atomnames = NTR_atomnames | set(env.atomnames[conf])

        confnames = [x.strip() for x in env.tpl[("CONFLIST", "CTR")].split(",")]
        CTR_atomnames = set()
        for conf in confnames:
            CTR_atomnames = CTR_atomnames | set(env.atomnames[conf])

        chains = []
        # count chains
        for res in self.residues:
            if res.chainid not in chains:
                chains.append(res.chainid)

        for chainid in chains:
            # get all residues of this chain
            aminoacids_in_chain = []
            others_in_chain = []
            for res in self.residues:
                if res.chainid == chainid:
                    if res.resname in AMINO_ACIDS:
                        aminoacids_in_chain.append(res)
                    else:
                        others_in_chain.append(res)

            if aminoacids_in_chain:
                ntr = aminoacids_in_chain[0]
                ctr = aminoacids_in_chain[0]
            else:
                continue

            for res in aminoacids_in_chain[1:]:
                if res.seqnum < ntr.seqnum:
                    ntr = res.resid[0]
                elif res.seqnum > ctr.seqnum:
                    ctr = res

            # verify bond for NTR
            atom_N =None
            for atom in ntr.conformers[0].atoms:
                if atom.name == " N  ":
                    atom_N = atom
                    break

            if not atom_N:
                logging.critical("No N atom found for NTR residue")
                return False

            ntr.flag = "ntr"
            for res in others_in_chain:
                if not ntr.flag:
                    break
                for conf in res.conformers:
                    if not ntr.flag:
                        break
                    for atom2 in conf.atoms:
                        d2 = d2vv(atom2.xyz, atom_N.xyz)
                        if d2 < clash_distance2:
                            ntr.flag = ""
                            break

            # verify bond for CTR
            atom_C =None
            for atom in ntr.conformers[0].atoms:
                if atom.name == " C  ":
                    atom_C = atom
                    break

            if not atom_C:
                logging.critical("No C atom found for CTR residue")
                return False

            ctr.flag = "ctr"
            for res in others_in_chain:
                if not ctr.flag:
                    break
                for conf in res.conformers:
                    if not ctr.flag:
                        break
                    for atom2 in conf.atoms:
                        d2 = d2vv(atom2.xyz, atom_C.xyz)
                        if d2 < clash_distance2:
                            ctr.flag = ""
                            break

        new_atoms = []
        for res in self.residues:
            for conf in res.conformers:
                for atom in conf.atoms:
                    if res.flag == "ntr" and atom.name in NTR_atomnames:
                        atom.resname = "NTR"
                        atom.resid = ("NTR", res.resid[1], res.resid[2], res.resid[3])
                    elif res.flag == "ctr" and atom.name in CTR_atomnames:
                        atom.resname = "CTR"
                        atom.resid = ("CTR", res.resid[1], res.resid[2], res.resid[3])
                    new_atoms.append(atom)

        self.load_atoms_single(new_atoms)

        return True


    def remove_exposed(self):


        return

    def sas(self):


        return

env = Env()
