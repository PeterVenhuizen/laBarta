#!/usr/bin/env python

"""
Create a spliced Salmon/KALLISTO transcriptome fasta file from a GTF file. 
"""

from natsort import natsorted
from fileParser import yield_fasta
import argparse
import subprocess

def gtf2fasta(gtf, genome_fasta):

	# Parse gtf and create exon bed
	with open('tmp.bed', 'w') as fout:
		for line in gtf:
			if not line.startswith('#'):
				c, source, feature, s, e, score, strand, phase, attributes = line.rstrip().split('\t')
				attr = {}
				for a in attributes.split(';'):
					if len(a):
						attr_name, attr_value = a.split(' "')
						attr[attr_name.strip()] = attr_value.replace('\"', '')

				if feature == "exon":
					fout.write( '{0}\t{1}\t{2}\t{3}_{0}:{1}-{2}_{4}\t1000\t{4}\n'.format(c, int(s)-1, e, attr['transcript_id'], strand) )

	# Get exon fasta
	cmd = "bedtools getfasta -name -s -fi {} -bed tmp.bed -fo tmp.fa".format(genome_fasta)
	subprocess.call(cmd, shell=True)

	# Parse exon fasta to spliced fasta
	trs, exon_dict = {}, {}
	for seq_record in yield_fasta('tmp.fa'):
		t_id, locus, other = seq_record.id.split('_')
		s, e = locus.split(':')[1].split('-')
		strand = other.split(':')[0]

		try: trs[t_id]['exons'].append((s, e))
		except KeyError: trs[t_id] = { 'exons': [(s, e)], 'strand': strand }

		exon_dict[s+'-'+e] = seq_record.seq

	# Output spliced transcriptome
	w = 70
	for t in natsorted(trs):

		ekeys = natsorted(trs[t]['exons']) if trs[t]['strand'] == '+' else natsorted(trs[t]['exons'], reverse=True)
		seq = ''.join([ exon_dict[s+'-'+e] for s, e in ekeys ])

		print '>{}'.format(t)
		for i in xrange(0, len(seq), w): print seq[i:i+w]

	subprocess.call('rm tmp.fa tmp.bed', shell=True)

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g', '--gtf', type=argparse.FileType('r'), default='-', help="Input GTF file")
	parser.add_argument('-f', '--genome-fasta', required=True, help="Genomic sequence fasta")
	args = parser.parse_args()

	gtf2fasta(args.gtf, args.genome_fasta)