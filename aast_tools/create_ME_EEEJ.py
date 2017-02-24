#!/usr/bin/env python

"""
Get all the exon-exon-junctions (EEJ), intron-exon-junctions (IEJ), exon-intron-junctions (EIJ),
and intron midpoints sequences based on the sets of transcript intron and exon fasta sequences. 
Calculate the number of unique mappable positions for each set of sequences. 
"""

from natsort import natsorted
import subprocess
import argparse
import gzip
import utils
import os
import re
import itertools

def get_bin(s, e, bin_size=100000):
	
	bin1 = int(s/bin_size)
	bin2 = int(e/bin_size)
	bin_code = str(bin1 * bin_size) + '-' + str((bin1+1) * bin_size)
	if bin1 == bin2: return [bin_code]
	else: return [ bin_code, str(bin2 * bin_size) + '-' + str((bin2+1) * bin_size) ]

def create_ME_EEEJ(exon_files, intron_files, output_dir, k=50, m=1, genome_FASTA="/home/venhuip8/data/master_transcriptome/a_thaliana/a_thaliana_TAIR10.fa"):

	# Parse exons and store all exons in bins of 5kb for filtering out
	# microexons overlapping with exons of the same or other genes. 
	exons, bins = {}, {}
	bin_size = 5000
	lociRE = re.compile('([0-5]):(\d+)-(\d+)')
	for f in exon_files:
		exons = utils.parse_fasta_file(exons, f)

		# Stuff exons in bins
		for line in open(f):
			if line.startswith('>'):
				try: 
					c, s, e = lociRE.search(line).groups()
					if c not in bins: bins[c] = {}
					s, e = natsorted(map(int, [s, e]))
					bin_codes = get_bin(s, e, bin_size)
					for bin_code in bin_codes:
						try: bins[c][bin_code].add((s, e))
						except KeyError: bins[c][bin_code] = set([(s, e)])
				except TypeError: pass

	# Parse introns
	introns = {}
	for f in intron_files:
		introns = utils.parse_fasta_file(introns, f)

	# Look for microexons and create EEEJs
	if not os.path.exists(output_dir): os.makedirs(output_dir)
	meRE = re.compile('(?=AG(.{3,15})GT)')
	ME_lib, ME_kmers = open(output_dir+'ME_lib.fa', 'w'), gzip.open(output_dir+'ME_lib.kmers.fa.gz', 'wt')
	for g_id in natsorted(introns):
	
		# Get orientation
		assert exons[g_id]['isFwd'] == introns[g_id]['isFwd'], "Exon and intron direction should be the same, quitting..."
		isFwd = introns[g_id].pop('isFwd')
		del exons[g_id]['isFwd']
	
		EEEJ = {} # Store in dictionary to avoid duplicates
		for t in natsorted(introns[g_id]):
	
			# Order introns and exons
			ikeys = natsorted(introns[g_id][t].keys()) if isFwd else natsorted(introns[g_id][t].keys(), reverse=True)
			ekeys = natsorted(exons[g_id][t].keys()) if isFwd else natsorted(exons[g_id][t].keys(), reverse=True)

			for i in xrange(len(ikeys)):
				seq = introns[g_id][t][ikeys[i]]
				s, e = ikeys[i].split('-')

				C1 = exons[g_id][t][ekeys[i]]
				C2 = exons[g_id][t][ekeys[i+1]]
				
				C1_e = ekeys[i].split('-')[1] if isFwd else ekeys[i].split('-')[0]
				C2_s = ekeys[i+1].split('-')[0] if isFwd else ekeys[i+1].split('-')[1]

				# Search for microexons ranging from 3 to 15nt. 
				# Ignore the original intron splice sites.
				# Microexons can overlap.
				isME = False
				for match in meRE.finditer(seq[2:-2]):
				
					me_seq = match.group(1)
					me_len = len(me_seq)
					me_start = int(s)+match.start()+3 if isFwd else int(e)-match.start()-4
					me_end = me_start+me_len if isFwd else me_start-me_len
					if not isFwd: me_start, me_end = me_end, me_start

					# Check if microexon is not within an annotated exon
					bin_codes = get_bin(me_start, me_end, bin_size)

					# Check all exons in bin
					keep = True
					for bin_code in bin_codes:
						if bin_code in bins[c] and any([ utils.is_overlapping(me_start, me_end, x, y) for (x, y) in bins[c][bin_code] ]): keep = False

					if keep: 

						isME = True
					
						exon_flank = 42-me_len
						ID = '>{}_{}_{}-{}_{}.INC'.format(g_id, C1_e, me_start, me_end, C2_s)
						EEEJ[ ID ] = '{}{}{}'.format(C1[-exon_flank:], me_seq, C2[:exon_flank])
		
				# Get the microexon exclusion EEJ
				if isME:
					
					EEEJ[ '>{}_{}-{}.EXC'.format(g_id, C1_e, C2_s) ] = C1[-42:] + C2[:42]

		for e in natsorted(EEEJ):
			ME_lib.write( '{}\n{}\n'.format(e, EEEJ[e]) )
		
			seq = EEEJ[e]
			N = (len(seq)-k)+1
			ME_kmers.write( ''.join([ '>{}|{}-{}\n{}\n'.format(g_id, i+1, i+k+1, seq[i:(i+k)]) for i in xrange(0, N, m) ]) )
			
	ME_lib.close(), ME_kmers.close()
	
	with open('tmp.sh', 'w') as fout:
		fout.write('#!/bin/bash\n')
		
		# Create TAIR10 + EEJ + EEEJ index
		cmd = 'bowtie-build {0},{1}ME_lib.fa {1}TAIR10_ME_lib &> /dev/null'.format(GENOME_FASTA, output_dir)
		fout.write(cmd+'\n')
		
		# Index EEJ + EEEJ library
		cmd = 'bowtie-build {0}ME_lib.fa {0}ME_lib &> /dev/null'.format(output_dir)
		fout.write(cmd+'\n')
		
		# Map to TAIR10 + EEJ + EEEJ
		cmd = """bowtie -f -p 6 -m 1 -v 2 --un {0}ME_lib.kmers.un.fa --max {0}ME_lib.kmers.max.fa {0}TAIR10_ME_lib <(zcat {0}ME_lib.kmers.fa.gz) | awk '{{ print ">"$1"\\n"$5 }}' > {0}ME_lib.kmers.mapped.fa""".format(output_dir)
		fout.write(cmd+'\n')
		
		# Map unique aligned kmers to EEJ + EEEJ and count effective number of unique mappable positions
		cmd = """bowtie -f -p 6 -m 1 -v 2 {0}ME_lib {0}ME_lib.kmers.mapped.fa | awk '{{ c[$3]++ }} END {{ OFS="\\t"; for (k in c) print k, c[k] }}' | sort > {0}ME_lib.eff""".format(output_dir)
		fout.write(cmd+'\n')
	
	subprocess.call('chmod +x tmp.sh', shell=True)
	subprocess.call('./tmp.sh', shell=True)
	subprocess.call('rm tmp.sh', shell=True)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-e', '--exon-fasta', nargs='+', required=True, help="Exon fasta file(s).")
	parser.add_argument('-i', '--intron-fasta', nargs='+', required=True, help="Intron fasta files(s).")
	parser.add_argument('-o', '--output-dir', required=True, help="Output directory.")
	parser.add_argument('-k', '--kmer', default=50, type=int, help="K-mer size (default = 50)")	
	parser.add_argument('-m', '--move-by', default=1, type=int, help="Move the sliding window by m nt (default = 1)")
	parser.add_argument('--GENOME-FASTA', default="/home/venhuip8/data/master_transcriptome/a_thaliana/a_thaliana_TAIR10.fa", help="Genome fasta file.")
	args = parser.parse_args()

	create_ME_EEEJ(args.exon_fasta, args.intron_fasta, utils.add_slash(args.output_dir), args.kmer, args.move_by, args.GENOME_FASTA)