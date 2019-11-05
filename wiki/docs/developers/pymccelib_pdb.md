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

### Read a PDB pdb file

### Read a MCCE pdb file

### Write a pdb file