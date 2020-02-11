# Protein structure manipulation

!!! Info
    After importing the library, you have access to the constants, attributes and functions:<br>
    ```>>> from pymccelib import *  ```

This part of functions handle the protein structure manipulation.

---
## Separate N and C terminals

Protein.identify_nc()

*Identify N and C terminus and put corresponding atoms as new residues.*

**Example:**
```
prot = Protein()
...
prot.identify_nc()
```

**Description:**

This subroutine does:

  * Starting from an existing protein structure, identify N and C terminus, take the backbone atoms that belong to N 
  or C terminus into a new residue.

The criteria for determining a N or C residue are *

  * a residue is an amino acid
  * sequence number is lowest or highest among the amino acids in a chain
  * For N terminus, N atom is not within bond distance of any atoms from other residues of the chain
  * For C terminus, O atom is not within bond distance of any atoms from other residues of the chain

This part of functions handle the protein structure manipulation.

---
## Solven Accessible Surface

Protein.sas()

*Calculate solvent accessible exposure of protein residue and conformer.*

**Example:**
```
prot = Protein()
...
prot.sas()
```

**Description:**

This subroutine does:

  * Calculate solvent exposed surface area, both absolute and fraction values, of each conformer.
  * Calculate solvent exposed surface area, both absolute and fraction values, of each residue.

Notes:

  * a surface point is considered accessible when the  point extended by the probe radius is not buried 
  * fraction sas is exposed area divided by exposed conformer or residue 
  * when calculating sas, other residues are assumed to have the native (1st) conformer only.



