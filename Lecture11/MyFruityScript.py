#!/usr/bin/python3
# Script to compare fruit
# Version 1, 25th Oct 2023
import os
import sys
# Allow for two fruit
# NOTE_incomplete: should probably include some error traps here for wrong input etc
first_fruit  = sys.argv[1]
second_fruit = sys.argv[2]
# Standard unchanging answer output
# NOTE_incomplete: a lot more could be done here...!
print(first_fruit.upper() + " are just as nice as  " + second_fruit.upper())
