# PDB IO

!!! Info
    After importing the library, you have access to the constants, attributes and functions:<br>
    ```>>> from pymccelib import *  ```

This part of functions handle the reading and writing of PDB files.

## MCCE protein structure organization
blockdiag {
P [label = "Protein"];
R1 [label = "residue 1"];
R2 [label = "residue 2"];
Rs [label = "residue ..."];
Rn [label = "residue N"];

C0 [label = "conformer backbone"];
C1 [label = "conformer 1"];
C2 [label = "conformer 2"];
Cs [label = "conformer ..."];
Cn [label = "conformer M"];

A1 [label = "atom 1"];
A2 [label = "atom 2"];
As [label = "atom ..."];
An [label = "atom n"];

A01 [label = "atom 1"];
A02 [label = "atom 2"];
A0s [label = "atom ..."];
A0n [label = "atom m"];

dots1 [shape = dots];
dots2 [shape = dots];
dots3 [shape = dots];
dots4 [shape = dots];
dots5 [shape = dots];
dots6 [shape = dots];
dots7 [shape = dots];

P -> R1;
P -> R2;
P -> Rs;
P -> Rn;

R1 -> C0;
R1 -> C1;
R1 -> C2;
R1 -> Cs;
R1 -> Cn;

R2 -> dots1;
R2 -> dots2;
Rs -> dots3;
Rn -> dots4;

C0 -> A01;
C0 -> A02;
C0 -> A0s;
C0 -> A0n;
C1 -> A1;
C1 -> A2;
C1 -> As;
C1 -> An;
C2 -> dots5;
Cs -> dots6;
Cn -> dots7;
}

The structure is organized as a hierarchical data structure as Protein -> Residue -> Conformer -> Atom.

Each residue contains one or multiple conformers. 

 * The first conformer is special. It contains atoms that are considered to be backbone, 
 which do not change charge and position in MCCE sampling.
 * Conformers other than the first conformer are candidates to compose a microstate, 
 which consists one conformer from each residue.

## Read/write a pdb
### MCCE PDB format
```
ATOM      8  CB  VAL A0003_001  17.786  38.593  41.279   2.000       0.000      01____M000
```

Column (0 based) | Definition | Format | Example 
--- | --- | --- | ---:
0-5 | Record name | %6s | ATOM 
6-10 | Serial number | %5d | 2 
12-15 | Atom name | %4s | CB 
17-19 | Residue name | %3s | VAL 
21 | Chain ID | %c | A
22-25 | Sequence number | %04d | 0003<sup>1</sup>
26 | Insertion code | %c | _<sup>2</sup>
27-29 | Conformer number | %03d | 001<sup>1</sup>
30-37 | X coordinate | %8.3f | 17.786
38-45 | Y coordinate | %8.3f | 38.593
46-53 | Z coordinate | %8.3f | 41.279
54-61 | PB radius | %8.3f | 2.000
66-73 | Charge | %8.3f | 0.000
80- | Conformer history | %s | 01____M000<sup>3</sup>

<sup>1</sup>: Filled with 0, and will later be used to compose conformer ID string.

<sup>2</sup>: If no insertion code, MCCE uses "_", and will later be used to compose conformer ID string.

<sup>3</sup>: Conformer making history:
  
  * Column 1-2: conformer type as a two letter code in residue CONFLIST in ftpl file. "BK", "01", "02", "+1" etc
  * Column 3: Rotamer making mechanism, a one letter code.
  * Column 4-6: Rotamer number, reindexed for each mechanism indicated by column 3.
  * Column 7: H atom placing mechanism, a one letter code.
  * Column 8-10: Rotamer number, reindexed for each mechanism indicated by column 7.

Rotamer making mechanism code explained, for column 3:
 
  * "O" = original 
  * "R" = rotated, sampled by repacking with neighboring residues
  * "X" = most exposed 
  * "H" = hydrogen-bond directed

Rotamer making mechanism code explained, for column 7:
 
  * "M" = torsion minimum
  * "H" = hydrogen-bond directed

---
### Read a native pdb file
Protein.load_nativepdb(pdb)

*Load native pdb file into Protein data structure.*

**Example:**
```
from pymccelib import *
 
env.init()

# get pdb file name from run.prm
pdbfile = env.prm("INPDB")

# initialize a protein object
prot = Protein()

# load pdb file
prot.load_nativepdb(pdbfile)
```

**Description:**

This subroutine does:

  * Read pdb file in native format (PDB format) and store in mcce protein hierarchical structure.

---
### Read a MCCE pdb file

---
### Write a pdb file
*Write mcce protein data structure to array of lines ready to output to stdout or file.*

**Example:**
```
prot=Protein()
...
lines = prot.pdblines()
sys.stdout.writelines(lines)
```

**Description:**

This subroutine does:

  * Read pdb file in native format (PDB format) and store in mcce protein hierarchical structure.
