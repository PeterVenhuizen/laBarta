#!/usr/bin/env python

"""
Create a spliced Salmon/KALLISTO transcriptome fasta file from a GTF file. 
"""

from natsort import natsorted
from fileParser import yield_fasta, parse_GTF
import argparse
import subprocess

def gtf2fasta(gtf_file, genome_fasta, select_feature="exon", t_id_attr="transcript_id", attr_sep=' "'):

	# Parse gtf and create exon bed
	gtf = parse_GTF(gtf_file, select_feature, t_id_attr, attr_sep, False)
	with open('tmp.bed', 'w') as fout:
		for t_id in natsorted(gtf):
			for s, e in gtf[t_id]['exons']:
				fout.write( '{0}\t{1}\t{2}\t{3}_{0}:{1}-{2}_{4}\t1000\t{4}\n'.format(gtf[t_id]['chr'], s-1, e, t_id, gtf[t_id]['strand']) )

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
	parser.add_argument('-g', '--gtf', required=True, help="Input GTF file")
	parser.add_argument('-f', '--genome-fasta', required=True, help="Genomic sequence fasta")
	parser.add_argument('--select-feature', default="exon", help="Feature to select in the GTF file")
	parser.add_argument('--transcript-attr', default='transcript_id', help="The transcript identifier attribute name")
	parser.add_argument('--attr-sep', default=' "', help="GFF attribute and attribute value separator")
	args = parser.parse_args()

	gtf2fasta(args.gtf, args.genome_fasta, args.select_feature, args.transcript_attr, args.attr_sep)