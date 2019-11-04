# Program Structure

## Main program:

### pymcce.py

This program works like the c version main program. It integrates the four steps in one program and use run.prm to to
 control which steps are executed.

blockdiag {
A [label = "run.prm", shape = "note"];
B [label = "param/*.ftpl", shape = "note"];
C [label = "pymcce.py", shape = ellipse];
D [label = "premcce.py", shape = ellipse];
E [label = "mkconf.py", shape = ellipse];
F [label = "mkenergy.py", shape = ellipse];
G [label = "monte.py", shape = ellipse];

A -> C;
B -> C;
C -> D;
C -> E; 
C -> F; 
C -> G; 
}

## Individual steps
Individual steps can run as standalone programs or be called by main program.

These steps exchange data by files.

### premcce.py: convert to mcce pdb

blockdiag {
A [label = "run.prm", shape = "note"];
B [label = "param/*.ftpl", shape = "note"];
B1 [label = "pdb", shape = "note"];
C [label = "pymcce.py", shape = ellipse];
D [label = "step1_out.pdb", shape = "note"];
E [label = "head1.lst", shape = "note"];

A -> C;
B -> C;
B1 -> C;
C -> D;
C -> E; 
}

Premcce.py reads a pdb file, and converts it to mcce pdb. In this process, the structure will have terminal residues 
and certain cofactors renamed to appropriate names so they can have proper ionization states and conformers.
  
Solvent and salt molecules are stripped off according to the rules in run.prm.
  
Head1.lst is a instruction file for the next step to make side chain ionization and conformation conformers. The 
default instruction is based on the rules in run.prm, and users can customize this file for more controls.

### mkconf.py: make conformers

blockdiag {
A [label = "run.prm", shape = "note"];
B [label = "param/*.ftpl", shape = "note"];
B1 [label = "step1_out.pdb", shape = "note"];
B2 [label = "head1.lst", shape = "note"];
C [label = "mkconf.py", shape = ellipse];
D [label = "step2_out.pdb", shape = "note"];
E [label = "head2.lst", shape = "note"];

A -> C;
B -> C;
B1 -> C;
B2 -> C
C -> D;
C -> E; 
}

Ionization and conformation conformers are recorded in step2_out.pdb. Head2.lst is formational, 
it records the rotamer making history and statistics. 

### mkenergy.py: caculate energy table

blockdiag {
A [label = "run.prm", shape = "note"];
B [label = "param/*.ftpl", shape = "note"];
B1 [label = "step2_out.pdb", shape = "note"];
C [label = "mkenergy.py", shape = ellipse];
D [label = "energies/*.opp", shape = "note"];
E [label = "head3.lst", shape = "note"];

A -> C;
B -> C;
B1 -> C;
C -> D;
C -> E; 
}

Mkenergy.py is the energy calculation step. It reads in step2_out.pdb and compute the side chain conformer pairwise 
interaction in energies/ directory, and side chain conformer self energy terms in head3.lst. 
    
    
### monte.py: Monte Carlo sampling of states

blockdiag {
A [label = "run.prm", shape = "note"];
B [label = "energies/*.opp", shape = "note"];
B1 [label = "head3.lst", shape = "note"];
C [label = "monte.py", shape = ellipse];
D [label = "fort.38", shape = "note"];
E [label = "mc_out", shape = "note"];
E [label = "microstates/*.ms", shape = "note"];

A -> C;
B -> C;
B1 -> C;
C -> D;
C -> E; 
}

Monte.py samples the side chain conformer occupancy therefore reports the probabilty of ionization at different 
conditions. Fort.38 is a conformer occupancy file, mc_out records the Monte Carlo progress including energy tracing, 
and microstates/*.ms are the microstates records in the sampling.

## Supporting modules:

### pymccelib.py: mcce modules
MCCE data structures and modules
    
### geometry.py: geometry operations
Geometry operations

## Analysis tools:
A collection fof analysis tools.

### fitpka.py: Titration curve fitting

### mfe.py: use mean field to analyze free energy

