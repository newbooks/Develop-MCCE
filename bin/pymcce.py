#!/usr/bin/env python

import logging

from pymccelib import *


def welcome():
    print("===========================================================")
    print("<<< MCCE Multi-Conformation Continuum Electrostatics >>>   ")
    print(" Marilyn Gunner's Lab at City College of New York, 2005    ")
    print("-----------------------------------------------------------")
    print("Version:        %s" % env.version)
    print("MCCE Home Page: https://sites.google.com/site/mccewiki     ")
    print("Support:        mgunner@ccny.cuny.edu                      ")
    print("Developed by:   Gunner Group @ City College of New York    ")
    print("Reference MCCE: If you publish data calculated with MCCE,  ")
    print("                please cite papers suggested in MCCE       ")
    print("                Home Page.                                 ")
    print("===========================================================")
    print("Last Updates:                                              ")
    print("   09/01/2019: Use free format tpl.                        ")
    print("===========================================================")
    print()
    return


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, format='%(message)s')
    welcome()
    env.init()

