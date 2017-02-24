#!/usr/bin/env python

"""
Identify events with a delta PSI/PIR/PSU above 15/25 across samples, if the sample std is not bigger than 5.
"""

import numpy as np
import argparse
import re
from natsort import natsorted

def parse_AS(f):

	d = {}
	for line in open(f):
		if not line.startswith('#') and 'Bn' not in line:
			gene_id, as_type, coord, psi, quality_scores, event_size = line.rstrip().split('\t')
			d[gene_id+'_'+coord+'_'+as_type] = { 'PSI': float(psi.split(':')[0]), 'event_size': event_size }
	return d

def QC(std):

	if std <= 1: score = "SOK"
	elif std <= 5: score = "OK"
	elif std <= 10: score = "B1"
	elif std <= 15: score = "B2"
	else: score = "Bn"

	return score

def calc_dPSI(s1_files, s2_files, dPSI=15, HQdPSI=25, max_std=5):

	# Parse AS files
	s1 = [ parse_AS(f) for f in s1_files ]
	s2 = [ parse_AS(f) for f in s2_files ]

	# Get shared events
	s1_shared = set.intersection(*[ set(s1[i].keys()) for i in xrange(len(s1)) ])
	s2_shared = set.intersection(*[ set(s2[i].keys()) for i in xrange(len(s2)) ])
	sample_intersection = s1_shared.intersection(s2_shared)

	print 'GENE\tCOORD\tEVENT_SIZE\tAS_TYPE\t{}\t{}\tDIRECTION\tDELTA'.format( '\t'.join([ f.split('/')[-1].split('.')[0] for f in s1_files ]), '\t'.join([ f.split('/')[-1].split('.')[0] for f in s2_files ]) )
	for event in natsorted(sample_intersection):

		# Get values
		s1_val, s2_val = [ [ s[i][event]['PSI'] for i in xrange(len(s)) ] for s in [s1, s2] ]

		# Get delta of means
		delta = abs(np.mean(s1_val)-np.mean(s2_val))

		# Check std and delta
		std = [ np.std(s1_val) <= max_std, np.std(s2_val) <= max_std ]

		if delta >= dPSI and all(std):
			g_id, coord, as_type = event.split('_')
			print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}'.format(g_id, coord, s[0][event]['event_size'], as_type, '\t'.join(map(str, s1_val)), '\t'.join(map(str, s2_val)), 'CTRL<TEST' if np.mean(s1_val) < np.mean(s2_val) else 'CTRL>TEST', delta)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-s1', required=True, help="Condition one files regex.", nargs='+')
	parser.add_argument('-s2', required=True, help="Condition two files regex.", nargs='+')
	parser.add_argument('--dPSI', default=10, type=int, help="Minimum delta PSI (default = 15).")
	parser.add_argument('--HQ-dPSI', default=25, type=int, help="High quality delta PSI (default = 25).")
	parser.add_argument('--max-std', default=5, type=int, help="Maximum internal condition variance.")
	args = parser.parse_args()

	calc_dPSI(args.s1, args.s2, args.dPSI, args.HQ_dPSI, args.max_std)