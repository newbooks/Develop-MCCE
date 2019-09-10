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
env.prm[TITR_PHD] = "1.0"
```
---

### Access env.prm

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
In run.prm, if a key is missing, or "DEFAULT" is used, the value will be interpreted in the context of the key:

Here is the default value list:

  * TPL_FOLDER
```
DEFAULT  tpl file folder path, "DEFAULT" to launch location (TPL_FOLDER)
```
&nbsp&nbsp guyg


---

## TPL file variables

## Methods in env
