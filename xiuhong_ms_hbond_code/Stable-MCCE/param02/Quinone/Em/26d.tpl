#2,6-dimethyl-Benzoquinone (Q-)        
#Agnes 7.9.08  
CONFLIST 26D        26DBK 26D01 26D-1 26DDM

NATOM    26DDM      0
NATOM    26DBK      0
NATOM    26D01      18
NATOM    26D-1      18

IATOM    26D01  C6  0
IATOM    26D01  C5  1
IATOM    26D01  H5  2
IATOM    26D01  H3  3
IATOM    26D01  C4  4
IATOM    26D01  C1  5
IATOM    26D01  C2  6
IATOM    26D01  C3  7
IATOM    26D01  CM6 8
IATOM    26D01  O1  9
IATOM    26D01  O4  10
IATOM    26D01  CM2 11
IATOM    26D01 2H1  12
IATOM    26D01 2H2  13
IATOM    26D01 2H3  14
IATOM    26D01 6H1  15
IATOM    26D01 6H2  16
IATOM    26D01 6H3  17
IATOM    26D-1  C6  0
IATOM    26D-1  C5  1
IATOM    26D-1  H5  2
IATOM    26D-1  H3  3
IATOM    26D-1  C4  4
IATOM    26D-1  C1  5
IATOM    26D-1  C2  6
IATOM    26D-1  C3  7
IATOM    26D-1  CM6 8
IATOM    26D-1  O1  9
IATOM    26D-1  O4  10
IATOM    26D-1  CM2 11
IATOM    26D-1 2H1  12
IATOM    26D-1 2H2  13
IATOM    26D-1 2H3  14
IATOM    26D-1 6H1  15
IATOM    26D-1 6H2  16
IATOM    26D-1 6H3  17



ATOMNAME 26D01    0  C6 
ATOMNAME 26D01    1  C5 
ATOMNAME 26D01    2  H5 
ATOMNAME 26D01    3  H3 
ATOMNAME 26D01    4  C4 
ATOMNAME 26D01    5  C1 
ATOMNAME 26D01    6  C2 
ATOMNAME 26D01    7  C3 
ATOMNAME 26D01    8  CM6 
ATOMNAME 26D01    9  O1 
ATOMNAME 26D01   10  O4 
ATOMNAME 26D01   11  CM2 
ATOMNAME 26D01   12 2H1
ATOMNAME 26D01   13 2H2
ATOMNAME 26D01   14 2H3
ATOMNAME 26D01   15 6H1
ATOMNAME 26D01   16 6H2
ATOMNAME 26D01   17 6H3
ATOMNAME 26D-1    0  C6
ATOMNAME 26D-1    1  C5
ATOMNAME 26D-1    2  H5
ATOMNAME 26D-1    3  H3
ATOMNAME 26D-1    4  C4
ATOMNAME 26D-1    5  C1
ATOMNAME 26D-1    6  C2
ATOMNAME 26D-1    7  C3
ATOMNAME 26D-1    8  CM6
ATOMNAME 26D-1    9  O1
ATOMNAME 26D-1   10  O4
ATOMNAME 26D-1   11  CM2
ATOMNAME 26D-1   12 2H1
ATOMNAME 26D-1   13 2H2
ATOMNAME 26D-1   14 2H3
ATOMNAME 26D-1   15 6H1
ATOMNAME 26D-1   16 6H2
ATOMNAME 26D-1   17 6H3




#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   26D01      0
PROTON   26D-1      0

PKA      26D01      0.0
PKA      26D-1      0.0

ELECTRON 26D01      0
ELECTRON 26D-1      1

EM       26D01      0.0
EM       26D-1      -80

RXN      26D01      -2.03
RXN      26D-1      -14.31 

#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  26D01  C6  sp2        0    C5   0    C1   0    CM6
CONNECT  26D01  C5  sp2        0    C6   0    C4   0    H5
CONNECT  26D01  H5  s          0    C5
CONNECT  26D01  H3  s          0    C3
CONNECT  26D01  C4  sp2        0    C5   0    C3   0    O4
CONNECT  26D01  C1  sp2        0    C6   0    C2   0    O1
CONNECT  26D01  C2  sp2        0    C1   0    C3   0    CM2
CONNECT  26D01  C3  sp2        0    C4   0    C2   0    H3
CONNECT  26D01  CM6 sp3        0    C6   0   6H1   0   6H2   0   6H3  
CONNECT  26D01  CM2 sp3        0    C2   0   2H1   0   2H2   0   2H3   
CONNECT  26D01  O1  s          0    C1
CONNECT  26D01  O4  s          0    C4
CONNECT  26D01 2H1  s          0    CM2
CONNECT  26D01 2H2  s          0    CM2
CONNECT  26D01 2H3  s          0    CM2
CONNECT  26D01 6H1  s          0    CM6
CONNECT  26D01 6H2  s          0    CM6
CONNECT  26D01 6H3  s          0    CM6

#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  26D-1  C6  sp2        0    C5   0    C1   0    CM6
CONNECT  26D-1  C5  sp2        0    C6   0    C4   0    H5
CONNECT  26D-1  H5  s          0    C5
CONNECT  26D-1  H3  s          0    C3
CONNECT  26D-1  C4  sp2        0    C5   0    C3   0    O4
CONNECT  26D-1  C1  sp2        0    C6   0    C2   0    O1
CONNECT  26D-1  C2  sp2        0    C1   0    C3   0    CM2
CONNECT  26D-1  C3  sp2        0    C4   0    C2   0    H3
CONNECT  26D-1  CM6 sp3        0    C6   0   6H1   0   6H2   0   6H3
CONNECT  26D-1  CM2 sp3        0    C2   0   2H1   0   2H2   0   2H3
CONNECT  26D-1  O1  s          0    C1
CONNECT  26D-1  O4  s          0    C4
CONNECT  26D-1 2H1  s          0    CM2
CONNECT  26D-1 2H2  s          0    CM2
CONNECT  26D-1 2H3  s          0    CM2
CONNECT  26D-1 6H1  s          0    CM6
CONNECT  26D-1 6H2  s          0    CM6
CONNECT  26D-1 6H3  s          0    CM6

#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   26D    C6  1.60
RADIUS   26D    C5  1.60
RADIUS   26D    H5  1.00
RADIUS   26D    H3  1.00
RADIUS   26D    C4  1.60
RADIUS   26D    C1  1.60
RADIUS   26D    C2  1.60
RADIUS   26D    C3  1.60
RADIUS   26D    O1  1.40
RADIUS   26D    O4  1.40
RADIUS   26D    CM6 1.60
RADIUS   26D    CM2 1.60
RADIUS   26D   2H1  1.00
RADIUS   26D   2H2  1.00
RADIUS   26D   2H3  1.00
RADIUS   26D   6H1  1.00
RADIUS   26D   6H2  1.00
RADIUS   26D   6H3  1.00


#NEUTRAL------
#23456789A123456789B123456789C
# opt ub3lyp/6-31++g(d,p)  nosymm pop=chelpg scf(maxcycle=200)   gk 7.31.06
CHARGE   26D01  O1  -0.46
CHARGE   26D01  C1   0.48
CHARGE   26D01  C4   0.71
CHARGE   26D01  O4  -0.54
CHARGE   26D01  C2   0.04
CHARGE   26D01  CM2 -0.05
CHARGE   26D01  C3  -0.38
CHARGE   26D01  CM6 -0.06
CHARGE   26D01  C5  -0.36
CHARGE   26D01  H3   0.17
CHARGE   26D01  C6   0.03
CHARGE   26D01  H5   0.17
CHARGE   26D01 2H1   0.04
CHARGE   26D01 2H2   0.04
CHARGE   26D01 2H3   0.04
CHARGE   26D01 6H1   0.04
CHARGE   26D01 6H2   0.04
CHARGE   26D01 6H3   0.04

# opt ub3lyp/6-31++g(d,p)  nosymm pop=chelpg scf(maxcycle=200)   Agnes 7.9.08
CHARGE   26D-1  O1  -0.65 
CHARGE   26D-1  C1   0.38
CHARGE   26D-1  C4   0.70
CHARGE   26D-1  O4  -0.75
CHARGE   26D-1  C2   0.02
CHARGE   26D-1  CM2  0.02
CHARGE   26D-1  C3  -0.47 
CHARGE   26D-1  CM6  0.01
CHARGE   26D-1  C5  -0.46
CHARGE   26D-1  H3   0.14 
CHARGE   26D-1  C6   0.02
CHARGE   26D-1  H5   0.13
CHARGE   26D-1 2H1  -0.01
CHARGE   26D-1 2H2  -0.01 
CHARGE   26D-1 2H3  -0.02
CHARGE   26D-1 6H1  -0.01
CHARGE   26D-1 6H2  -0.01 
CHARGE   26D-1 6H3  -0.02 
 
#TORSION  26D    HO2  O2   C2   C3   f      1.800         2    180.00

#ParaNam|Res  |Atom|Param/toggle
TRANS    26D          t

