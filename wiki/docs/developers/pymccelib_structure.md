# Protein structure manipulation

!!! Info
    After importing the library, you have access to the constants, attributes and functions:<br>
    ```>>> from pymccelib import *  ```

This part of functions handle the protein structure manilulation.

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


