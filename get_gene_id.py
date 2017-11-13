#!/usr/bin/env python

"""
Get the gene_id value from a line
"""

import re
import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-f', type=argparse.FileType('r'), default='-', help="Whatever file containing a gene_id \".*\"; value.")
args = parser.parse_args()

gene_id = re.compile('gene_id "(.*?)";')
for line in args.f:
	try: print re.search(gene_id, line).group(1)
	except AttributeError: pass