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

def score_seq(seq, PWM):
	return sum([ log((PWM[i][base]/0.25)+0.0001, 2) for i, base in enumerate(seq) if base in PWM[i].keys() ])

def wholeSeqPWM(background_fasta, soi_fasta, PWM_files, min_score=65):
	
	# Parse PWM files
	for f in PWM_files:
		PWM = { 'pwm': read_PWM_from_file(f), 'min': 0, 'max': 0 }
		k = len(PWM['pwm']) # Get PWM size
		
		kmers = {}
		# Iterate over sequences and score the kmers in the same
		# size as the PWMs.
		
		# Only use the background sequences for the PWM score scaling
		for record in yield_fasta(background_fasta):
			N = (len(record.seq)-k)+1
			for i in xrange(0, N, 1):
				
				score = score_seq(record.seq[i:i+k], PWM['pwm'])
				PWM['min'] = score if score < PWM['min'] else PWM['min']
				PWM['max'] = score if score > PWM['max'] else PWM['max']
		
		# Score and store the kmers in the SOIs and save the 
		# sequence identifiers matching the kmers.
		for record in yield_fasta(soi_fasta):
			
			N = (len(record.seq)-k)+1
			for i in xrange(0, N, 1):
				
				kmer = record.seq[i:i+k]
				
				try:
					kmers[kmer]['IDs'].add( record.id )
				except KeyError: 

					# Score and store
					score = score_seq(kmer, PWM['pwm'])
					kmers[kmer] = { 'score': score, 'IDs': set([ record.id ]) }

					# Update min and max scores
					PWM['min'] = score if score < PWM['min'] else PWM['min']
					PWM['max'] = score if score > PWM['max'] else PWM['max']
					
		# Rescale the scores and output hits above the min_score cutoff
		for kmer in natsorted(kmers):
			
			if kmer['score'] > 0: norm_score = ((50*kmer['score'])/PWM['max'])+50
			else: norm_score = -(50*(kmer['score']-PWM['min'])/PWM['min'])
				
			if norm_score >= min_score:
				pass

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('--soi', required=True, help="Fasta file containing the sequences of interest.")
	parser.add_argument('-b', '--background', required=True, help="Fasta file with 'background' sequence to compare against.")
	parser.add_argument('--PWM', required=True, nargs='+', help="PWM files (output from get_PWMs.py) or similarly structured files.")
	parser.add_argument('--min-score', type=int, default=65, help="Minimal PWM score cutoff for sequence reporting")
	args = parser.parse_args()
	
	
