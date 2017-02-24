#!/usr/bin/env python

"""
Get all the exon-exon junctions (EEJs) by combining all the donor and acceptor sites in the 
annotated transcript models. Calculate the number of unique mappable positions for each EEJ.
"""

from natsort import natsorted
import subprocess
import argparse
import gzip
import utils
import os

def create_SSB_EEJ(exon_files, output_dir, k=50, m=1, GENOME_FASTA="/home/venhuip8/data/master_transcriptome/a_thaliana/a_thaliana_TAIR10.fa"):

	# Parse exon fasta
	seqs = {}
	for f in exon_files:
		seqs = utils.parse_fasta_file(seqs, f)

	# Get EEJs and split to kmers
	if not os.path.exists(output_dir): os.makedirs(output_dir)
	EEJ_lib, EEJ_kmers = open('{}SSB_lib.fa'.format(output_dir), 'w'), gzip.open('{}SSB_lib.kmers.fa.gz'.format(output_dir), 'wt')
	for g_id in natsorted(seqs):
		
		EEJ = {} # Store in dict to avoid duplicate EEJs
		exons = seqs[g_id]
		isFwd = exons.pop('isFwd')
		
		for t in natsorted(exons):
			ekeys = natsorted(exons[t].keys()) if isFwd else natsorted(exons[t].keys(), reverse=True)
			for i in xrange(len(ekeys)):
				x = ekeys[i]
				e = x.split('-')[1] if isFwd else x.split('-')[0]
				for y in ekeys[i+1:]:
					s = y.split('-')[0] if isFwd else y.split('-')[1]
					EEJ[ '{}_{}-{}'.format(g_id, e, s) ] = exons[t][x][-42:] + exons[t][y][:42]
		
		# Write EEJs and kmers
		for eej in natsorted(EEJ):
			EEJ_lib.write( '>{}\n{}\n'.format(eej, EEJ[eej]) )
			
			seq = EEJ[eej]
			N = (len(seq)-k)+1
			EEJ_kmers.write( ''.join([ '>{}|{}-{}\n{}\n'.format(g_id, i+1, i+k+1, seq[i:(i+k)]) for i in xrange(0, N, m) ]) )
					
	EEJ_lib.close(), EEJ_kmers.close()
	
	with open('tmp.sh', 'w') as fout:
		fout.write('#!/bin/bash\n')
	
		# Create TAIR10 + EEJ index
		cmd = 'bowtie-build {0},{1}SSB_lib.fa {1}TAIR10_SSB_lib &> /dev/null'.format(GENOME_FASTA, output_dir)
		fout.write(cmd+'\n')

		# Index EEJ library
		cmd = 'bowtie-build {0}SSB_lib.fa {0}SSB_lib &> /dev/null'.format(output_dir)
		fout.write(cmd+'\n')
		
		# Map to TAIR10 + EEJ
		cmd = """bowtie -f -p 6 -m 1 -v 2 --un {0}SSB_lib.kmers.un.fa --max {0}SSB_lib.kmers.max.fa {0}TAIR10_SSB_lib <(zcat {0}SSB_lib.kmers.fa.gz) | awk '{{ print ">"$1"\\n"$5 }}' > {0}SSB_lib.kmers.mapped.fa""".format(output_dir)
		fout.write(cmd+'\n')
		
		# Map unique aligned kmers to EEJ and count effective number of unique mappable positions
		cmd = """bowtie -f -p 6 -m 1 -v 2 {0}SSB_lib {0}SSB_lib.kmers.mapped.fa | awk '{{ c[$3]++ }} END {{ OFS="\\t"; for (k in c) print k, c[k] }}' | sort > {0}SSB_lib.eff""".format(output_dir)
		fout.write(cmd+'\n')
	
	subprocess.call('chmod +x tmp.sh', shell=True)
	subprocess.call('./tmp.sh', shell=True)
	subprocess.call('rm tmp.sh', shell=True)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-e', '--exon-fasta', nargs='+', required=True, help="Exon fasta file(s).")
	parser.add_argument('-o', '--output-dir', required=True, help="Output directory.")
	parser.add_argument('-k', '--kmer', default=50, type=int, help="K-mer size (default = 50)")	
	parser.add_argument('-m', '--move-by', default=1, type=int, help="Move the sliding window by m nt (default = 1)")
	parser.add_argument('--GENOME-FASTA', default="/home/venhuip8/data/master_transcriptome/a_thaliana/a_thaliana_TAIR10.fa", help="Genome fasta file.")
	args = parser.parse_args()

	create_SSB_EEJ(args.exon_fasta, utils.add_slash(args.output_dir), args.kmer, args.move_by, args.GENOME_FASTA)