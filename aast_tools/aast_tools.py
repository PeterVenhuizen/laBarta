#!/usr/bin/env python

"""
Arabidopsis Alternative Splicing and Transcription Tools (AAST-TOOLS) 
a.k.a. VAST-TOOLS for Arabidopsis wrapper
"""

import os
import argparse
import subprocess
from utils import add_slash

def aast_tools(datasets, reads_dir, lib_dir, output_dir="./aast_tools/", scripts_dir="/home/venhuip8/scripts/aast_tools/"):

	tmp_dir = "{}tmp/".format(output_dir)
	count_dir = "{}module_count/".format(output_dir)
	psi_dir = "{}PSI/".format(output_dir)
	log_file = "{}aast_tools.log".format(output_dir)

	for f in [tmp_dir, count_dir, psi_dir, psi_dir+"SS/", psi_dir+"ES_TB/", psi_dir+"ES_SSB/", psi_dir+"IR/", psi_dir+"ME_v2/"]:
		if not os.path.exists(f): os.makedirs(f)

	for d in ["TB", "SSB", "ME_v2", "IR_IEJ", "IR_MID"]: 
		if not os.path.exists(count_dir+d): os.makedirs(count_dir+d)

	for line in open(datasets):
		condition, fq1, fq2 = line.rstrip().split('\t')

		# Create bash script
		with open("{}tmp.sh".format(output_dir), 'w') as fout:

			print(condition)
			fout.write( '#!/bin/bash\n' )
			fout.write( 'echo {} >> {}\n'.format(condition, log_file) )

			fout.write( '# Split RNAseq to 50nt kmers\n' )
			fout.write( 'python {}fastq2kmers.py -i {}{} -o {}fwd.fa.gz\n'.format(scripts_dir, reads_dir, fq1, tmp_dir) )
			fout.write( 'python {}fastq2kmers.py -i {}{} -o {}rev.fa.gz\n'.format(scripts_dir, reads_dir, fq2, tmp_dir) )

			fout.write( '# Map to TAIR10\n' )
			fout.write( '(bowtie -f -p 6 -m 1 -v 2 --un {0}unmapped.fa --max {0}max.fa {1}TAIR10 <(zcat {0}fwd.fa.gz),<(zcat {0}rev.fa.gz) | awk \'{{ print ">"$1"\\n"$6 }}\' > {0}TAIR10_mapped.fa) 2> {2}\n'.format(tmp_dir, lib_dir, log_file) )

			fout.write( '# Transcript based module: Map unmapped to EEJ\n' )
			fout.write( '(bowtie -f -p 6 -m 1 -v 2 --un {0}un.fa --max {0}max.fa {1}transcript_based/TB_lib {0}unmapped.fa | awk \'{{ c[$4]++ }} END {{ OFS="\\t"; for (k in c) print k, c[k] }}\' | sort > {2}TB/{3}.txt) 2> {4}\n'.format(tmp_dir, lib_dir, count_dir, condition, log_file) )

			fout.write( '# Splice site based module: Map unmapped to EEJ\n' )
			fout.write( '(bowtie -f -p 6 -m 1 -v 2 --un {0}un.fa --max {0}max.fa {1}splice_site_based/SSB_lib {0}unmapped.fa | awk \'{{ c[$4]++ }} END {{ OFS="\\t"; for (k in c) print k, c[k] }}\' | sort > {2}SSB/{3}.txt) 2> {4}\n'.format(tmp_dir, lib_dir, count_dir, condition, log_file) )

			fout.write( '# Microexon module: Map unmapped to EEEJ (-m 1 -v 2)\n' )
			fout.write( '(bowtie -f -p 6 -m 1 -v 2 --un {0}un.fa --max {0}max.fa {1}microexon/ME_lib {0}unmapped.fa | awk \'{{ c[$4]++ }} END {{ OFS="\\t"; for (k in c) print k, c[k] }}\' | sort > {2}ME_v2/{3}.txt) 2> {4}\n'.format(tmp_dir, lib_dir, count_dir, condition, log_file) )

			fout.write( '# Intron retention module: Map against TAIR10 + EEJ\n' )
			fout.write( '(bowtie -f -p 6 -m 1 -v 2 --un {0}un.fa --max {0}max.fa {1}intron_retention/TAIR10_IR.EEJ_lib <(zcat {0}fwd.fa.gz),<(zcat {0}rev.fa.gz) | awk \'{{ print ">"$1"\\n"$6 }}\' > {0}mapped.fa) 2> {2}\n'.format(tmp_dir, lib_dir, log_file) )
			fout.write( '# Intron retention module: Map unique to EEJ + IEJ\n' )
			fout.write( '(bowtie -f -p 6 -m 1 -v 2 --un {0}un.fa --max {0}max.fa {1}intron_retention/IR.IEJ_lib {0}mapped.fa | awk \'{{ c[$3]++ }} END {{ OFS="\\t"; for (k in c) print k, c[k] }}\' | sort > {2}IR_IEJ/{3}.txt) 2> {4}\n'.format(tmp_dir, lib_dir, count_dir, condition, log_file) )
			fout.write( '# Intron retention module: Map TAIR10 unique to midpoint sequences\n' )
			fout.write( '(bowtie -f -p 6 -m 1 -v 2 --un {0}un.fa --max {0}max.fa {1}intron_retention/IR.MID_lib {0}TAIR10_mapped.fa | awk \'{{ c[$3]++ }} END {{ OFS="\\t"; for (k in c) print k, c[k] }}\' | sort > {2}IR_MID/{3}.txt) 2> {4}\n'.format(tmp_dir, lib_dir, count_dir, condition, log_file) )

			fout.write( '(python {0}calculate_PSI.py -f {1}AtRTD2.exons.fa -e {1}splice_site_based/SSB_lib.eff -c {2}SSB/{3}.txt -AS SS > {4}SS/{3}.txt) 2> {4}SS/{3}_cov.txt\n'.format(scripts_dir, lib_dir, count_dir, condition, psi_dir) )
			fout.write( 'python {0}calculate_PSI.py -f {1}AtRTD2.exons.fa -e {1}splice_site_based/SSB_lib.eff -c {2}SSB/{3}.txt -AS ES > {4}ES_SSB/{3}.txt\n'.format(scripts_dir, lib_dir, count_dir, condition, psi_dir) )
			fout.write( 'python {0}calculate_PSI.py -f {1}AtRTD2.exons.fa -e {1}transcript_based/TB_lib.eff -c {2}TB/{3}.txt -AS ES > {4}ES_TB/{3}.txt\n'.format(scripts_dir, lib_dir, count_dir, condition, psi_dir) )
			fout.write( '(python {0}calculate_PSI.py -e {1}microexon/ME_lib.eff -c {2}ME_v2/{3}.txt -AS ME > {4}ME_v2/{3}.txt) 2> {4}ME_v2/{3}_cov.txt\n'.format(scripts_dir, lib_dir, count_dir, condition, psi_dir) )
			fout.write( '(COUNT=$(mktemp -q /tmp/COUNT.XXXXXX); cat {0}IR_*/{1}.txt > $COUNT; EFF=$(mktemp -q /tmp/EFF.XXXXXX); cat {2}intron_retention/IR*.eff > $EFF; python {3}calculate_PSI.py -e $EFF -c $COUNT -AS IR > {4}IR/{1}.txt) 2> {4}IR/{1}_cov.txt\n'.format(count_dir, condition, lib_dir, scripts_dir, psi_dir) )

		subprocess.call( 'chmod +x {0}tmp.sh; {0}./tmp.sh'.format(output_dir), shell=True )

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-d', '--datasets', required=True, help="Text file containing datasets information.")
	parser.add_argument('-r', '--reads-dir', required=True, help="Reads directory.")
	parser.add_argument('-l', '--lib', default="/home/venhuip8/scripts/aast_tools/lib/AtRTD2/", help="AAST-TOOLS exon-exon junction libraries")
	parser.add_argument('-o', '--output-dir', default="./aast_tools/", help="AAST-TOOLS output directory.")
	parser.add_argument('-s', '--scripts', default="/home/venhuip8/scripts/aast_tools/", help="AAST-TOOLS scripts location.")
	args = parser.parse_args()

	aast_tools(args.datasets, args.reads_dir, args.lib, add_slash(args.output_dir), add_slash(args.scripts))