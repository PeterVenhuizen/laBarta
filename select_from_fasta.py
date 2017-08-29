#!/usr/bin/env python

"""
	Select from fasta. Extract the fasta records which (exactly)
	match the identifiers given in a separate file.
"""

from __future__ import print_function
import sys
import argparse
from fileParser import yield_fasta

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

def select_from_fasta(fasta_file, selection_file, modus="exact"):

	selects = [ line.rstrip() for line in open(selection_file) ]

	for record in yield_fasta(fasta_file):
		if modus == "exact":
			if record.id in selects: print('>{}\n{}'.format(record.id, record.seq))
			else: eprint('>{}\n{}'.format(record.id, record.seq))
		elif modus == "contains":
			looking = True
			for s in selects:
				if s in record.id: 
					print('>{}\n{}'.format(record.id, record.seq))
					looking = False
					break

			if looking: 
				eprint('>{}\n{}'.format(record.id, record.seq))

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('--fasta', required=True, help="Fasta file to extract records from.")
	parser.add_argument('--selection', required=True, help="Selection file containing identifiers.")
	parser.add_argument('--modus', choices=["exact", "contains"], default="exact", help="Extract records which exactly match the identifiers or when the identifier is contained in the fasta record identifier.")
	args = parser.parse_args()

	select_from_fasta(args.fasta, args.selection, args.modus)