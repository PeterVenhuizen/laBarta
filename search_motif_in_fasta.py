#!/usr/bin/env python

"""
Search a motif in the input fasta file. 
"""

from __future__ import division
from fileParser import yield_fasta
from sequence import Sequence
import argparse
import re

def search_motifs(fasta, motif, name, reverse=False):

	name = name.replace(' ', '_')

	# Convert motif to regex if necessary
	if any([ 
		'B' in motif, 'D' in motif, 'H' in motif,
		'I' in motif, 'K' in motif, 'M' in motif,
		'N' in motif, 'R' in motif, 'S' in motif,
		'V' in motif, 'W' in motif, 'Y' in motif
	]): motif = Sequence(motif).iupac2regex()
	p = re.compile('(?=(%s))' % (motif)) # Allows for overlapping matches

	total_len, fwd_hits, rev_hits, out = 0, {}, {}, []
	for record in yield_fasta(fasta):
		seq_len = len(record.seq)
		total_len += seq_len
		
		# Forward search
		for m in p.finditer(record.seq):
			out.append( '{}\t{}\t{}\t{}|{}\t1000\t+\n'.format(record.id, m.start(), m.start()+len(m.group(1)), name, m.group(1)) )
			try: fwd_hits[m.group(1)] += 1
			except KeyError: fwd_hits[m.group(1)] = 1

		# Reverse search
		if reverse: 
			seq = Sequence(record.seq).get_reverse_complement()
			for m in p.finditer(seq):
				out.append( '{}\t{}\t{}\t{}|{}\t1000\t-\n'.format(record.id, seq_len-(m.start()+len(m.group(1))), seq_len-m.start(), name, m.group(1)) )
				try: rev_hits[m.group(1)] += 1
				except KeyError: rev_hits[m.group(1)] = 1
	
	# Summarize results
	fwd_total = sum([ fwd_hits[h] for h in fwd_hits ])
	rev_total = sum([ rev_hits[h] for h in rev_hits ])

	header = '#Search space: {}nt\n#Motif: {} ({})\n'.format(total_len, motif, name)
	header += '#FORWARD\n#Total\t{}\t{:.3f}/kb\n'.format( fwd_total, fwd_total * (1000/total_len) )
	for h in fwd_hits: header += '#{}\t{}\t{:.3f}/kb\n'.format( h, fwd_hits[h], fwd_hits[h] * (1000/total_len) )

	if reverse:
		header += '#REVERSE\n#Total\t{}\t{:.3f}/kb\n'.format( rev_total, rev_total * (1000/total_len) )
		for h in rev_hits: header += '#{}\t{}\t{:.3f}/kb\n'.format( h, rev_hits[h], rev_hits[h] * (1000/total_len) )

	print header, ''.join(out)
	
if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-f', '--fasta', help="Sequences in fasta format", required=True)
	parser.add_argument('-m', '--motif', help="Motif sequence in IUPAC or regex format", required=True)
	parser.add_argument('--name', help="Motif name", required=False, nargs='+', default='')
	parser.add_argument('--reverse', dest='reverse', action='store_true', default=False)
	args = parser.parse_args()
	
	search_motifs(args.fasta, args.motif, '_'.join(args.name), args.reverse)