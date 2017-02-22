#!/usr/bin/env python

'''
	Name: get_AS_landscape.py
	Author: Petrus Andreas Venhuizen
	E-mail: peter.venhuizen@univie.ac.at
	Institute: Max F. Perutz Laboratories / Medical University of Vienna
	Description: Report the alternative splicing landscape of a set of transcripts.
	Required packages: natsort (https://pypi.python.org/pypi/natsort)
	
	Created on: 01-10-2015
	Last changed: 22-02-2017
	
	To-do:
	- 	Change parse_gtf(). Allow for missing exon_number attribute.
	-	Add gff file support
'''

from __future__ import print_function
import os
import argparse
from natsort import natsorted, versorted
from fileParser import parse_GTF

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

def go_landscaping(transcripts, ref_id, chromosome, event_type, detailed=False):
	pass

def run(gtf_file, ref, event_type, output_dir):
	
	if not os.path.exists(output_dir): os.makedirs(output_dir)

	# Parse transcripts in the gtf_file

	# Parse reference transcripts

if __name__ == '__main__':

	AS_choices = ['ES', 'SS', 'IR', 'FL']

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g', '--gtf', required=True, help="Transcripts in GTF format.")
	parser.add_argument('-r', '--reference-trs', required=True, help="List of reference transcripts for the input GTF genes.")
	parser.add_argument('-e', '--event-type', nargs='+', choices=AS_choices, default=AS_choices, help="AS events to look for.")
	parser.add_argument('-o', '--output', default="./", help="Output path/file")
	parser.add_argument('--version', action='version', version='v0.3.0')
	args = parser.parse_args()