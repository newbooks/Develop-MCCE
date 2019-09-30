# PyMCCE Library

!!! Info
    After importing the library, you have access to the constants, attributes and functions:<br>
    ```>>> from pymccelib import *  ```
    
## Constants
These are global variables

---
### ROOMT

**Definition:** 

Room temperature.

**Example:**
```
>>> from pymccelib import *
>>> print(ROOMT)
298.15
```
---    
### PH2KCAL

**Definition:**

Conversion from 1 pH unit to kcal/mol.
 
**Example:**
```
>>> from pymccelib import *
>>> print(PH2KCAL)
1.364
```
--- 

### KCAL2KT
**Definition:**

Conversion from 1 kcal/mol to KT

**Example:**
```
>>> from pymccelib import *
>>> print(KCAL2KT)
1.688
```
---
### KJ2KCAL
**Definition:**

Conversion from 1kj/mol to kcal/mol.

**Example:**
```
>>> from pymccelib import *
>>> print(KJ2KCAL)
0.239
```

---

## MCCE run environment
MCCE run environment can be retrieved from dictionary env.prm. Both key and values are in string type. 

---
### run.prm --> env.prm
The file run.prm is translated to env.prm line by line. The last string in parenthesis is the key, 
and the first string is the value.

In file run.prm, line
```
f        step 1: pdb -> mcce pdb                            (DO_PREMCCE)
```

translates to:
```
env.prm["DO_MCCE"] = "f"
```

and line
```
1.0      pH interval                                        (TITR_PHD)
```
translates to:
```
env.prm["TITR_PHD"] = "1.0"
```
---

### Accessing env.prm

**Example:**

To access env.prm as a dictionary.

```
>>> from pymccelib import *     # import 
>>> env.init()                  # load run.prm and ftpl files
>>> print(env.prm)              # print out env.prm
{'INPDB': 'prot.pdb', 'DO_PREMCCE': 'f', 'DO_ROTAMERS': 'f', 'DO_ENERGY': 'f', 'DO_MONTE': 'f', 'TPL_FOLDER': '/home/jmao/projects/Develop-MCCE/param', 'EXTRA': './extra.ftpl', 'TITR_TYPE': 'ph', 'TITR_PH0': '0.0', 'TITR_PHD': '1.0', 'TITR_EH0': '0.0', 'TITR_EHD': '30.0', 'TITR_STEPS': '15', 'BIG_PAIRWISE': '5.0', 'MONTE_FLIPS': '3', 'MONTE_T': '298.15', 'MONTE_NITER': '2000', 'MONTE_RUNS': '6', 'NSTATE_MAX': '1000000'}
```

To access individual parameter:
```
>>> print(env.prm["DO_PREMCCE"])
f
```

---
### Exceptions and defaults 

*In run.prm, if a key is missing, or "DEFAULT" is used, the value will be interpreted in the context of the key.*


**Example:**
```
DEFAULT  tpl file folder path, "DEFAULT" to launch location (TPL_FOLDER)
```

The key TPL_FOLDER will be assigned a default value.


---
#### TPL_FOLDER
  
Topology file folder. Default is relative location "```../param```" to the executable location.

blockdiag {
A0 [label = "mcce root/"];
A[label = "bin/"];
B[label = "pymcce.py"];
C[label = "param/"];
D[label = "*.ftpl"];

A0 -> A -> B;
A0 -> C -> D;
}


The above example demonstrates where pymcce.py looks for topology files *.ftpl. 

--- 

#### DELPHI_EXE
  
The path of PB solver delphi. Default is the same directory where the main executable resides. 

blockdiag {
A0 [label = "mcce root/"];
A[label = "bin/"];
B[label = "pymcce.py"];
C[label = "delphi"];

A0 -> A -> B;
A -> C;
}


The above example demonstrates where delphi executable is if pymcce.py was the main program. 

--- 
#### Scaling factors 
Scaling factors are used by step 4, microstate sampling.

|Name | Definition | Default|
|---|---|---|
| SCALING_VDW0 | Scaling factor of intra-conformer VDW | "1.0" |
| SCALING_VDW1 | Scaling factor of sidechain to backbone VDW | "1.0" |
| SCALING_VDW | Scaling factor of sidechain to sidechain VDW | " 1.0" |
| SCALING_TORS | Scaling factor of torsion energy | " 1.0" |
| SCALING_ELE | Scaling factor of electrostatic interaction | " 1.0" |
| SCALING_DSOLV | Scaling factor of desolvation energy | " 1.0" |

---


## Topology records

These records hold definitions of ligands, conformers, and molecular parameters. The records are located as 
dictionary in env.tpl.

### ftpl -> env.tpl

Topology file is now free format. It has up to three keys, and a string as value. Key and value are 
separated by ":", and three parts of key are separated by ",".

**Example:**

```
CONNECT, " N  ", GLUBK: sp2, " ?  ", " CA ", " H  "
```

After converted, it is:
```
env.tpl[(CONNECT, " N  ", GLUBK)] = 'sp2, " ?  ", " CA ", " H  "'
```

Rules of reading ftpl files:

  * **Comments:** Anything after "#" will not be read.
  * **Free format:** Spaces are ignored unless inside quotes.
  * **Delimiter of key and value:** ":" divides key and value.
  * **Delimiter of key parts:** "," divides key parts. There can be up to 3 parts of a key.
  * **Enclosing space in key and value:** space can be key or value when enclosed in quotes. 
  * **Records with identical key:** the later key-value pair will overwrite the previous, with a warning.
  * **Order of reading ftpl files:** ```ftpl files in param/``` folder are read in alphabetic order and ```extra
  .ftpl``` defined in run.prm is read last. 
  
### Native and derived tpl entries
Native entries are directly translated from ftpl files. All values are strings.

#### Native entries

| Primary Key | Key example | Explanation |
|---|---|---|
| CONFLIST | ("CONFLIST", "GLU") | returns a conformer list from a residue name |
| CONNECT |("CONNECT", " N  ", "GLUBK")| returns orbital type and a list of connected atoms |
| CHARGE |("CHARGE", "GLUBK", " N  ") | returns atom charge |
| RADIUS |("RADIUS", "GLUBK",  " N  ")| dielectric boundary radius, VDW radius and energy well depth|
| CONFORMER|("CONFORMER", "GLU01")| Conformer self energy terms |
| ROTATE |("ROTATE", "GLU")| Residue rotatable bonds |

!!! Note
    Native entries are stored in env.tpl{}.

#### Derived entries
Derived entries are converted or structured entries derived from native entries. 

| Name | Example | Explanation |
| --- | --- | --- |
|env.atomenames{}| env.atomnames["GLUBK"] | atom names derived from CONNECT records|
|env.confnames{} | env.confnames["GLU"] | conformer names of a residue |

!!! Note
    Derived entries are stored in separate env attributes.

### Accessing env.tpl
Native and derived tpl entries can be accessed by visiting respected dictionaries in env. Value data type varies.

**Example**
```
>>> from pymccelib import *
>>> env.init()
>>> print(env.tpl[("CHARGE", "GLUBK", " N  ")])
-0.350
>>> print(env.atomnames["GLU01"])
[' CB ', ' HB2', ' HB3', ' CG ', ' HG2', ' HG3', ' CD ', ' OE1', ' HE1', ' OE2']
```


## Functions in env
---
### env.init()
*Initialize mcce running environment.*

**Synopsis:**
```
env.init()
```

**Description:**

This subroutine does:

  * Read run.prm, and set default values for missing parameters.
  * Read param/*.ftpl and extra.ftpl; initialize env.tpl{} and derived tpl entries.  

---

### env.print_runprm()
*Print all env.prm entries*

**Synopsis:**

```
env.print_runprm()
```

This subroutine prints out all env.prm entries for debug purpose.
 
---

### env.print_scaling()
*Print all scaling factors*

**Synopsis:**

```
env.print_scaling()
```

This subroutine prints out scaling factors used by microstate sampling for debug purpose.
 
---

### load_ftpl(self, file)
*Load an individual ftpl file*

**Synopsis:**

```
env.load_ftpl(ftpl_file_path)
```

This subroutine loads an individual ftpl file. MCCE loads all ftple files from location defined by line (TPL_FOLDER) in run.prm. One can load additional ftpl files to supplement or overwrite entries previously defined. 
 
---

### ftpl2tpl(self, file)
*Convert env.tpl (in ftpl format) to a tpl file for backward compatibility*

**Synopsis:**

```
env.ftpl2tpl(tpl_file_name)
```

This subroutine converts existing env.tpl{} to a single tpl file that old mcce can read.
---
