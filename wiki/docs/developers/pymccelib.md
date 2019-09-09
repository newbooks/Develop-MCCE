# PyMCCE Library

!!! Info
    After importing the library, you have access to the constants, variables and functions:<br>
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

### Exceptions and defaults 

## TPL file variables
