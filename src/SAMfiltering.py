#!/usr/bin/env python

import sys
import os, errno
import argparse
import time
import pysam

start_time = time.time()

parser = argparse.ArgumentParser(description='Filtering catalog from contigs alignments. Script was writen by C. Pouchon (2021).')
parser.add_argument("-p","--pop", help="population file",
                    type=str)
parser.add_argument("-i","--input", help="alignment file in sam format",
                    type=str)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
