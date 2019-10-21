#! /usr/bin/python3

import os
import glob

def main():
    tests = []

    for otherpath in glob.iglob(os.path.join(os.getcwd(), "*.txt")):
        with open(otherpath, "r") as rna_file:
            individuals = rna_file.readlines()[4:]
            tests += individuals
            print(otherpath + " was read.")

    with open("all_individuals.txt", "x") as all_individuals:
        for line in tests:
            all_individuals.write(line)

if __name__ == "__main__":
    main()
