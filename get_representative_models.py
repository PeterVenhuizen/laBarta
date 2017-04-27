#!/usr/bin/env python

"""
	Select for each the gene the transcript with the longest protein
	coding isoform. If there are multiple transcripts of the same 
	length for a particular gene, select the one with the highest 
	expression in WT (SR protein wild-type control) as determined by 
	Salmon.
"""

from natsort import natsorted
from collections import Counter
from fileParser import yield_fasta, parse_salmon
import argparse
import operator

def get_representative_models(protein_fasta, salmon_files):

	# Get the lengths of all isoforms per gene
	fa = {}
	for record in yield_fasta(protein_fasta):
		g_id = record.id.split('.')[0]
		try: fa[g_id][record.id] = len(record.seq)
		except KeyError: fa[g_id] = { record.id: len(record.seq) }

	# Read in Salmon files
	salmon = [ parse_salmon(f) for f in salmon_files ]

	# Select the longest isoform for each gene
	for g_id in natsorted(fa):

		longest = max(fa[g_id].iteritems(), key=operator.itemgetter(1))[1]
		t_ids = [ t_id for t_id in fa[g_id] if fa[g_id][t_id] == longest ]
		
		# If only one transcript is found, output directly
		if len(t_ids) == 1: print ''.join(t_ids)
		
		# Select the transcript with the highest expression,
		# based on majority voting
		else: 
			highest = Counter()
			for s in salmon:
				expr = { t_id: s[g_id][t_id] for t_id in t_ids }
				highest[max(expr.iteritems(), key=operator.itemgetter(1))[0]] += 1
			print highest.most_common(1)[0][0]
			

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-p', '--protein', required=True, help="Fasta file with isoform protein translations.")
	parser.add_argument('-s', '--salmon', nargs='+', help="Salmon TPM files.")
	args = parser.parse_args()

	get_representative_models(args.protein, args.salmon)