#!/usr/bin/env python

"""
Whole sequence PWM scoring
"""

from __future__ import division
from natsort import natsorted
from fileParser import yield_fasta
from math import log
import argparse

def read_PWM_from_file(f):
	
	pwm = []
	parse = True
	bases = []
	for line in open(f):
		
		if not line.startswith('#'):
			
			# Only parse the PWM counts and calculate the fractions,
			# this way we avoid any potential rounding inaccuracies
			if line in ['\n', '\r\n']: parse = False
			elif parse:
				
				cols = line.rstrip().split('\t')
				try: 
					total = sum([ float(x) for x in cols[1:] ])
					pwm.append({ bases[i-1]: float(cols[i])/total for i in xrange(1, len(cols)) })
				except ValueError:
					bases = cols[1:]
	
	return pwm

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-f', '--fasta', required=True, help="Fasta input file")
	parser.add_argument('--PWM', required=True, nargs='+', help="PWM files (output from get_PWMs.py) or similarly structured files.")
	parser.add_argument('--min-score', type=int, default=65, help="Minimal PWM score cutoff for sequence reporting")
	args = parser.parse_args()
	
	
