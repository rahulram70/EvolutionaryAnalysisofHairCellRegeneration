#!/usr/bin/env python3

import os
import sys

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    res_dir = os.path.abspath(os.path.join(script_dir, os.pardir))

    # go up 2 directories to project dir
    for i in range(2):
        res_dir = os.path.abspath(os.path.join(res_dir, os.pardir))

    # Setup paths to needed directories
    #
    pl_fl = script_dir + "/gather_domain_info.pl"
    res_dir += "/resources/data-raw/protein-sequences/"

    #  add all tr - pr associations to a list
    for file in os.listdir(res_dir):
        if (".txt" in file):
            prg = "perl {}".format(pl_fl)
            command = "{} '{}' &".format(prg, file)
            print(command)
            os.system(command)

if __name__ == '__main__':
    main()
