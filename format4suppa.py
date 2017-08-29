#!/usr/bin/env python

"""
format4suppa.py -> Format the Salmon/Kallisto data into the 
expression file format required by SUPPA.
"""

import argparse
from natsort import natsorted
from fileParser import yield_salmon, yield_kallisto

def format4suppa(files, sample_names, file_type="salmon"):
	
	if file_type not in ["salmon", "kallisto"]:
		print "ERROR: Unrecognized file type. Quitting..."
	else:

		if file_type == "salmon":
			expr = [ { record['transcript_id']: record['tpm'] for record in yield_salmon(f) } for f in files ]
		elif file_type == "kallisto":
			expr = [ { record['transcript_id']: record['tpm'] for record in yield_kallisto(f) } for f in files ]

		print '\t'.join(sample_names)
		for t_id in natsorted(expr[0]):
			print '{}\t{}'.format(t_id, '\t'.join([ str(x[t_id]) for x in expr ]))

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-f', '--files', nargs='+', required=True, help="File(s) to reformat.")
	parser.add_argument('-n', '--names', nargs='+', required=True, help="Samples names.")
	parser.add_argument('-t', '--type', choices=["salmon", "kallisto"], default="salmon", help="File format.")
	args = parser.parse_args()

	format4suppa(args.files, args.names, args.type)