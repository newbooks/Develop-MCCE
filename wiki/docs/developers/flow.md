# Python MCCE

A Python version of Multi-Conformation Continuum Electrostatics simulation program to model protein electrostatic interactions.

## Program Mechanism
A protein or protein complex structure is broken down into pieces, usually amino acid residues and cofactors. Each piece may have charge, protonation, and position conformers. These peices interact with each other, and the piece conformers are sampled by Monte Carlo sampling to find their probabilty of occupancy at various conditions, such as pH, redox potention, chemical potential, and introduced charges.

## Program Structure

### Step 1: Break structue into residues

``` mermaid
graph TD
A1(run.prm)-->B1((read prm))
A2(input.pdb)-->C1((read pdb))
A3(tpl files)-->B4((read tpl))
subgraph Step 0: Load run control and tpl files
    B1-->B2[run_prm db]
    B1-->B4
    B4-->B5[tpl db]
end
B2-->C1
B5-->C1
subgraph Step 1: Break into residues
    C1-->C2((verify pdb))
    C2-->C3((to residues))
    C3-->C4((missing atoms))
    C4-->C5((ligand bonds))
    C5-->C6((NTR, CTR))
end
C6 -->C7(step1_out.pdb)
C6 -->C8(head1.lst)
```

### Step 2: Make residue conformers

``` mermaid
graph TD
A1(run.prm)-->B1((read prm))
A3(tpl files)-->B4((read tpl))
subgraph Step 0: Load run control and tpl files
    B1-->B2[run_prm db]
    B1-->B4
    B4-->B5[tpl db]
end
B2-->C1
B5-->C1
B6(step1_out.pdb)-->C1
B7(head1.lst)-->C1((readin))
subgraph Step 2: Make residue conformers
    C1-->C2((charge<br/>conformers))
    C2-->C3((position<br/> conformers))
    C3-->C4((hbond <br/>conformers))
end
C4 -->C5(step2_out.pdb)
C4 -->C6(head2.lst)
```

### Step 3. Compute energy table
``` mermaid
graph TD
A1(run.prm)-->B1((read prm))
A3(tpl files)-->B4((read tpl))
subgraph Step 0: Load run control and tpl files
    B1-->B2[run_prm db]
    B1-->B4
    B4-->B5[tpl db]
end
B2-->C1
B6(step2_out.pdb)-->C1((readin))
subgraph Step 3: Make residue conformers
	C1-->C2((RXN delphi))
	C2-->C3((PW delphi))
	C3-->C4((boundary<br/>correction))
end
C4-->C5(energies/*.opp)
C4-->C6(head3.lst)
```

### Step 4. Monte Carlo sampling: accessible states
```  mermaid
graph TB
A1(energies/*.opp)-->B1(readin)
A2(head3.lst)-->B1
A3(run.prm)-->A4((read prm))
A4-->A5[prm db]
A5-->B1
B1-->B2((group into residues))
B2-->B3((verify flags))
B3-->B4((big list))
B4-->B5((loop over pH))
B5-->B6((independent<br\>runs))
B6-->B7((Monte Carlo sampling))
B7-->B7.1((record state))
B7.1-->B8{over?}
B8--no-->B6
B8--yes-->B9{next pH?}
B9--no-->B5
B9--yes-->B10((end))
```

