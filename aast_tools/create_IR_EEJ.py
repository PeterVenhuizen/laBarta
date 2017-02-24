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

def split_to_kmer(ID, seq, k=50, m=1):

	N = (len(seq)-k)+1
	return ''.join([ '>{}|{}-{}\n{}\n'.format(ID, i+1, i+k+1, seq[i:(i+k)]) for i in xrange(0, N, m) ])

def create_IR_EEJ(exon_files, intron_files, output_dir, k=50, m=1, GENOME_FASTA="/home/venhuip8/data/master_transcriptome/a_thaliana/a_thaliana_TAIR10.fa"):

	# Parse fasta
	exons = {}
	for f in exon_files:
		exons = utils.parse_fasta_file(exons, f)

	introns = {}
	for f in intron_files:
		introns = utils.parse_fasta_file(introns, f)

	# Get EEJ, IEJ and midpoint-sequences
	if not os.path.exists(output_dir): os.makedirs(output_dir)
	fEEJ = open( '{}IR.EEJ_lib.fa'.format(output_dir), 'w' )
	fIEJ = open( '{}IR.IEJ_lib.fa'.format(output_dir), 'w' )
	fMID = open( '{}IR.MID_lib.fa'.format(output_dir), 'w' )
	EEJ_kmers, MID_kmers = gzip.open( '{}IR.EEJ_lib.kmers.fa.gz'.format(output_dir), 'wt' ), gzip.open( '{}IR.MID_lib.kmers.fa.gz'.format(output_dir), 'wt' )
	for g_id in natsorted(exons):

		# Check if gene has introns
		if g_id in introns: 

			# Get orientation
			assert exons[g_id]['isFwd'] == introns[g_id]['isFwd'], "Exon and intron direction should be the same, quitting..."
			isFwd = introns[g_id].pop('isFwd')
			del exons[g_id]['isFwd']

			IEJ, EEJ, MID = {}, {}, {}
			for t in introns[g_id]:

				# Order introns and exons
				ikeys = natsorted(introns[g_id][t].keys()) if isFwd else natsorted(introns[g_id][t].keys(), reverse=True)
				ekeys = natsorted(exons[g_id][t].keys()) if isFwd else natsorted(exons[g_id][t].keys(), reverse=True)
			
				for i, x in enumerate(ikeys):
				
					E1 = exons[g_id][t][ekeys[i]]
					I1 = introns[g_id][t][x]
					E2 = exons[g_id][t][ekeys[i+1]]
				
					s_old = ekeys[i].split('-')[1]
					e_old = ekeys[i+1].split('-')[0]
					s = ekeys[i].split('-')[1] if isFwd else ekeys[i].split('-')[0]
					e = ekeys[i+1].split('-')[0] if isFwd else ekeys[i+1].split('-')[1]

					# Get exon-exon junction
					EEJ[ '{}_EEJ_{}-{}'.format(g_id, s, e) ] = E1[-42:] + E2[:42]

					# Get exon-intron and intron-exon junctions
					IEJ[ '{}_EIJ_{}-{}'.format(g_id, s, e) ] = E1[-42:] + I1[:42]
					IEJ[ '{}_IEJ_{}-{}'.format(g_id, s, e) ] = I1[-42:] + E2[:42]

					# Get intron midpoints
					in_len = len(I1)
					if in_len <= 200: 
						MID[ '{}_MID_{}-{}'.format(g_id, s, e) ] = I1
					else:
						mid = int(round(in_len/2))
						MID[ '{}_MID_{}-{}'.format(g_id, s, e) ] = I1[(mid-100):(mid+100)]

			# Output EEJs
			fIEJ.write( ''.join([ '>{}\n{}\n'.format(s, IEJ[s]) for s in IEJ ]) )
			fEEJ.write( ''.join([ '>{}\n{}\n'.format(s, EEJ[s]) for s in EEJ ]) )
			fMID.write( ''.join([ '>{}\n{}\n'.format(s, MID[s]) for s in MID ]) )

			# Get kmers
			for eej in IEJ: EEJ_kmers.write( split_to_kmer(g_id, IEJ[eej], k, m) )
			for eej in EEJ: EEJ_kmers.write( split_to_kmer(g_id, EEJ[eej], k, m) )
			for mid in MID: MID_kmers.write( split_to_kmer(g_id, MID[mid], k, m) )
	
	EEJ_kmers.close(), MID_kmers.close()
	fIEJ.close(), fEEJ.close(), fMID.close()
	
	# Get unique mappable positions of IEJ and EEJ
	with open('tmp.sh', 'w') as fout:
		fout.write('#!/bin/bash\n')
		
		# Create TAIR10 index
		fout.write( 'bowtie-build {} {}{}\n'.format(GENOME_FASTA, output_dir, '.'.join(GENOME_FASTA.split('/')[-1].split('.')[:-1])) )

		# Create TAIR10 + EEJ index
		fout.write( 'bowtie-build {0},{1}IR.EEJ_lib.fa {1}TAIR10_IR.EEJ_lib &> /dev/null\n'.format(GENOME_FASTA, output_dir) )

		# Index EEJ library
		fout.write( 'bowtie-build {0}IR.EEJ_lib.fa,{0}IR.IEJ_lib.fa {0}IR.IEJ_lib &> /dev/null\n'.format(output_dir) )
		
		# Map to TAIR10 + EEJ
		fout.write( """bowtie -f -p 6 -m 1 -v 2 --un {0}IR.IEJ_lib.kmers.un.fa --max {0}IR.IEJ_lib.kmers.max.fa {0}TAIR10_IR.EEJ_lib <(zcat {0}IR.EEJ_lib.kmers.fa.gz) | awk '{{ print ">"$1"\\n"$5 }}' > {0}IR.IEJ_lib.kmers.mapped.fa\n""".format(output_dir) )
		
		# Map unique aligned kmers to EEJ and count effective number of unique mappable positions
		fout.write( """bowtie -f -p 6 -m 1 -v 2 {0}IR.IEJ_lib {0}IR.IEJ_lib.kmers.mapped.fa | awk '{{ c[$3]++ }} END {{ OFS="\\t"; for (k in c) print k, c[k] }}' | sort > {0}IR.IEJ_lib.eff\n""".format(output_dir) )
	
	subprocess.call('chmod +x tmp.sh', shell=True)
	subprocess.call('./tmp.sh', shell=True)
	subprocess.call('rm tmp.sh', shell=True)
	
	# Get unique mappable positions of the intron midpoints
	with open('tmp.sh', 'w') as fout:
		fout.write('#!/bin/bash\n')
		
		# Map kmers to TAIR10
		fout.write( """bowtie -f -p 6 -m 1 -v 2 --un {0}IR.MID_lib.kmers.un.fa --max {0}IR.MID_lib.kmers.max.fa {0}{1} <(zcat {0}IR.MID_lib.kmers.fa.gz) | awk '{{ print ">"$1"\\n"$5 }}' > {0}IR.MID_lib.kmers.mapped.fa\n""".format(output_dir, '.'.join(GENOME_FASTA.split('/')[-1].split('.')[:-1])) )
	
		# Index midpoint library
		fout.write( 'bowtie-build {0}IR.MID_lib.fa {0}IR.MID_lib &> /dev/null\n'.format(output_dir) )
		
		# Map unique aligned kmers to midpoint library and count effective number of mappable positions
		fout.write( """bowtie -f -p 6 -m 1 -v 2 {0}IR.MID_lib {0}IR.MID_lib.kmers.mapped.fa | awk '{{ c[$3]++ }} END {{ OFS="\\t"; for (k in c) print k, c[k] }}' | sort > {0}IR.MID_lib.eff\n""".format(output_dir) )

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

	create_IR_EEJ(args.exon_fasta, args.intron_fasta, utils.add_slash(args.output_dir), args.kmer, args.move_by, args.GENOME_FASTA)