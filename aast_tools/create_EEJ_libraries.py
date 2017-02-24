#!/usr/bin/env python

"""
Wrapper for EEJ libraries creation. See individual EEJ library scripts for more details:
	* Transcript based => create_TB_EEJ.py
	* Splice site based => create_SSB_EEJ.py
	* Microexon => create_ME_EEEJ.py
	* Intron retention => create_IR_EEJ.py
"""

import os
import argparse
import subprocess
from utils import add_slash

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-e', '--exon-fasta', nargs='+', required=True, help="Exon fasta file(s).")
parser.add_argument('-i', '--intron-fasta', nargs='+', required=True, help="Intron fasta files(s).")
parser.add_argument('-o', '--output-dir', required=True, help="Output directory.")
parser.add_argument('-k', '--kmer', default=50, type=int, help="K-mer size (default = 50)")	
parser.add_argument('-m', '--move-by', default=1, type=int, help="Move the sliding window by m nt (default = 1)")
parser.add_argument('--GENOME-FASTA', default="/home/venhuip8/data/master_transcriptome/a_thaliana/a_thaliana_TAIR10.fa", help="Genome fasta file.")
args = parser.parse_args()

args.output_dir = add_slash(args.output_dir)

# TRANSCRIPT BASED EEJ
output_path = "{}transcript_based/".format(args.output_dir)
if not os.path.exists(output_path): os.makedirs(output_path)
cmd = "python create_TB_EEJ.py -e {} -o {} -k {} -m {} --GENOME-FASTA {}".format(args.exon_fasta, output_path, args.kmer, args.move_by, args.GENOME_FASTA)
subprocess.call(cmd, shell=True)

# SPLICE SITE BASED EEJ
output_path = "{}splice_site_based/".format(args.output_dir)
if not os.path.exists(output_path): os.makedirs(output_path)
cmd = "python create_SSB_EEJ.py -e {} -o {} -k {} -m {} --GENOME-FASTA {}".format(args.exon_fasta, output_path, args.kmer, args.move_by, args.GENOME_FASTA)
subprocess.call(cmd, shell=True)

# MICROEXON EEEJ
output_path = "{}microexon/".format(args.output_dir)
if not os.path.exists(output_path): os.makedirs(output_path)
cmd = "python create_ME_EEEJ.py -e {} -i {} -o {} -k {} -m {} --GENOME-FASTA {}".format(args.exon_fasta, args.intron_fasta, output_path, args.kmer, args.move_by, args.GENOME_FASTA)
subprocess.call(cmd, shell=True)

# INTRON RETENTION EEJ
output_path = "{}intron_retention/".format(args.output_dir)
if not os.path.exists(output_path): os.makedirs(output_path)
cmd = "python create_IR_EEJ.py -e {} -i {} -o {} -k {} -m {} --GENOME-FASTA {}".format(args.exon_fasta, args.intron_fasta, output_path, args.kmer, args.move_by, args.GENOME_FASTA)
subprocess.call(cmd, shell=True)