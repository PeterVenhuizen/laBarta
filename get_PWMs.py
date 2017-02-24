#!/usr/bin/env python

"""
Reassemble unspliced transcript sequences from the exon and intron 
fasta files and generate position weight matrices. Intron and exon
bed files can be obtained from GTF with gtf2exon_intron.py
"""

from __future__ import division
from natsort import natsorted
from fileParser import yield_fasta
import numpy as np
import subprocess
import argparse
import os

def parse_fasta_file(fasta_file):
	
	#>32899360|exon-1 | x-x | scaffold_350:2-1060 FORWARD LENGTH=1059::scaffold_350:1-1060(+)

	d = {}
	for record in yield_fasta(fasta_file):

		t_id = record.id.split('|')[0]
		locus = record.id.split(' ')[4].split(':')[1]
		isFwd = True if "FORWARD" in record.id.upper() else False

		try: d[t_id][locus] = record.seq
		except KeyError: d[t_id] = { locus: record.seq, 'isFwd': isFwd, 'chr': record.id.split(' ')[4].split(':')[0] }

	return d

def update_PWM(pwm, seq):

	for i, base in enumerate(seq):
		if base in ['A', 'C', 'G', 'T']:
			if pwm[i] == None:
				pwm[i] = { 'A': 0, 'C': 0, 'G': 0, 'T': 0 }
			pwm[i][base] += 1
	return pwm

def print_PWM(pwm, title, output_dir):

	bases = ['A', 'C', 'G', 'T']
	len_pwm = len(pwm)
	positions = xrange(-3, 10) if len_pwm == 13 else xrange(-14, 3)
	positions = [ x+1 if x >= 0 else x for x in positions ]

	# Create figure
	with open("save_bar_graph.R", 'w') as fout:
		fout.write( "#!/usr/bin/env Rscript\n" )
		fout.write( "library(ggplot2)\n" )
		fout.write( "library(scales)\n\n" )
	
		fout.write( "df <- data.frame(\n" )
		fout.write( "\tbases=c(\t{}),\n".format(', '.join([ '"' + x + '"' for x in bases * len_pwm ])) )
		fout.write( "\tposition=c({}),\n".format( ', '.join([ ', '.join([ '"'+str(y)+'"' for y in [x]*len(bases) ]) for x in positions ]) ) )
		fout.write( "\tpercent=c(\n" )

		for i in xrange(len_pwm):
			mean = sum([ pwm[i][b] for b in bases ])
			fout.write( '\t\t{}'.format( ', '.join([ str(pwm[i][b]/mean) for b in bases ]) ) )
			if i+1 != len_pwm: 
				fout.write( ',\n' )
			else: 
				fout.write( '\n\t)\n)' )

		fout.write( "\ndf$position <- factor(df$position, levels=c({}))".format(', '.join([ '"'+str(i)+'"' for i in positions ])) )

		fout.write( "\npng(file=\"{}/{}.png\", width={}, height=500)".format(output_dir, title, len_pwm*54) )
		fout.write( "\nggplot(data=df, aes(x=position, y=percent, fill=bases)) +" )		
		fout.write( "\n\tggtitle(\"{}\") + theme(plot.title = element_text(hjust = 0.5)) +".format(title) )
		fout.write( "\n\tgeom_bar(stat=\"identity\") +" )
		fout.write( "\n\tgeom_text(aes(label=ifelse(percent >= 0.05, levels(bases), \"\")), position=position_stack(vjust=0.5), color=\"white\") +" )
		#fout.write( "\n\tscale_fill_manual(values=c(\"#fef4ae\", \"#bac2f2\", \"#ffaca8\", \"#a8f3ba\")) +" )
		fout.write( "\n\tscale_fill_manual(values=c(\"forestgreen\", \"blue\", \"darkorange\", \"red\")) +" )
		fout.write( "\n\tlabs(x=\"\", y=\"\", fill=\"\") +")
		fout.write( "\n\tscale_y_continuous(labels=percent)" )

	subprocess.call("Rscript save_bar_graph.R", shell=True)
	subprocess.call("rm save_bar_graph.R", shell=True)

	# Output PWM text file
	with open('{}/{}.txt'.format(output_dir, title), 'w') as fout:
		fout.write( '\t' + '\t'.join(bases) + '\n')
		for i in xrange(len_pwm):
			fout.write( "{}\t{}\n".format(positions[i], '\t'.join([ str(pwm[i][b]) for b in bases ])) )

		fout.write( '\n\t' + '\t'.join(bases) + '\n')
		for i in xrange(len_pwm):
			total = sum([ pwm[i][b] for b in bases ])
			fout.write( "{}\t{}\n".format(positions[i], '\t'.join([ "{:.4f}".format(pwm[i][b]/total) for b in bases ])) )

def run(exon_fasta, intron_fasta, output_dir="./"):

	# Parse fasta
	exons = parse_fasta_file(exon_fasta)
	introns = parse_fasta_file(intron_fasta)

	PWMs = {}
	uniq_snowflakes = set()
	for t_id in natsorted(exons):

		if t_id in introns:

			# Check if gene has introns
			assert exons[t_id]['isFwd'] == introns[t_id]['isFwd'], "Exon and intron direction should be the same, quitting..."
			isFwd = introns[t_id].pop('isFwd')
			del exons[t_id]['isFwd']
			c = introns[t_id].pop('chr')
			del exons[t_id]['chr']

			# Order introns and exons
			ikeys = natsorted(introns[t_id].keys()) if isFwd else natsorted(introns[t_id].keys(), reverse=True)
			ekeys = natsorted(exons[t_id].keys()) if isFwd else natsorted(exons[t_id].keys(), reverse=True)

			for i, x in enumerate(ikeys):

				# Get upstream and downstream exons plus the intron
				E1 = exons[t_id][ekeys[i]]
				I = introns[t_id][x]
				E2 = exons[t_id][ekeys[i+1]]

				# Get the 5' and 3' splice sites
				five = I[:2]
				three = I[-2:]
				dinu = '{}-{}'.format(five, three)

				locus = '{}:{}'.format(c, x)
				if locus not in uniq_snowflakes:

					uniq_snowflakes.add(locus)
					five_seq = E1[-3:] + I[:10]
					three_seq = I[-14:] + E2[:3]
					if len(five_seq) == 13 and len(three_seq) == 17: 

						# Generate all PWM down below...
						if dinu not in PWMs:
							PWMs[dinu] = { '5SS': [None]*13, '3SS': [None]*17, 'N': 0 }
						PWMs[dinu] = { 
							'5SS': update_PWM(PWMs[dinu]['5SS'], five_seq),
							'3SS': update_PWM(PWMs[dinu]['3SS'], three_seq),
							'N': PWMs[dinu]['N'] + 1
						}

	if not os.path.exists(output_dir): os.makedirs(output_dir)
	for dinu in PWMs:
		if PWMs[dinu]['N'] >= 10:
			print_PWM(PWMs[dinu]['5SS'], '{}_{}'.format(dinu, '5SS'), output_dir)
			print_PWM(PWMs[dinu]['3SS'], '{}_{}'.format(dinu, '3SS'), output_dir)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-e', '--exon-fasta', help="Comma separated list of exon fasta files", required=True)
	parser.add_argument('-i', '--intron-fasta', help="Comma separated list of intron fasta files", required=True)
	parser.add_argument('-o', '--output-dir', default="./", help="Output directory")
	args = parser.parse_args()

	run(args.exon_fasta, args.intron_fasta, args.output_dir)