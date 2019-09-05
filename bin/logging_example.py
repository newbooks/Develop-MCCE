#!/usr/bin/env python
import logging


def welcome():
    logging.info("Welcome from function!")
    return


if __name__ == "__main__":
    # Set logging level. DEBUG, INFO, WARNING, ERROR and CRITICAL are in the order from more details to less
    # This has to be declared only once! and before any logging messages!
    logging.basicConfig(level=logging.DEBUG, format='%(message)s')
    logging.info("Welcome from main program!")
    welcome()
