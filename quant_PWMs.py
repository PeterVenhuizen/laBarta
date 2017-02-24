#!/usr/bin/env python

"""
Calculate the splice site strength for a set of sequences
and the species background of all introns. 
"""

from __future__ import division
from natsort import natsorted
from fileParser import yield_fasta
from math import log
import argparse

def parse_fasta_file(fasta_file):
	
	#>32899360|exon-1 | x-x | scaffold_350:2-1060 FORWARD LENGTH=1059::scaffold_350:1-1060(+)

	d = {}
	for seq_record in yield_fasta(fasta_file):

		t_id = seq_record.id.split('|')[0]
		locus = seq_record.id.split(' ')[4].split(':')[1]
		isFwd = True if "FORWARD" in seq_record.id.upper() else False

		try: d[t_id][locus] = seq_record.seq
		except KeyError: d[t_id] = { locus: seq_record.seq, 'isFwd': isFwd, 'chr': seq_record.id.split(' ')[4].split(':')[0] }

	return d

def read_PWM_from_file(f):
	
	pwm = []
	parse = True
	bases = []
	for line in open(f):

		if not line.startswith('#'):

			# Only parse the counts PWM and calculate the percentages
			# I think this is better, so we avoid rounding errors/inaccuracies 
			if line in ['\n', '\r\n']: parse = False 
			elif parse:

				cols = line.rstrip().split('\t')
				try:
					total = sum([ float(x) for x in cols[1:] ])
					pwm.append( { bases[i-1]: float(cols[i])/total for i in xrange(1, len(cols)) } )
				except ValueError:
					bases = cols[1:]

	return pwm

def score_splice_site(seq, PWM):

	return sum([ log((PWM[i][base]/0.25)+0.0001, 2) for i, base in enumerate(seq) if base in PWM[i].keys() ])

def update_min_max(seq, pwm):

	score = score_splice_site(seq, pwm['pwm'])
	pwm['min'] = score if score < pwm['min'] else pwm['min']
	pwm['max'] = score if score > pwm['max'] else pwm['max']
	return pwm

def rescale_score(PWM_score, ss_min, ss_max):
    ''' Scale the PWM LOD score to 0-100 
    
        ((b-a)*(PWM_score-min)/(max-min))+a
        
        If the PWM score is positive -> a = 50, b = 100, min = 0
        If the PWM score is negative -> a = 0, b = 50, max = 0
    '''

    if PWM_score > 0: norm_score = ((50*PWM_score)/ss_max)+50
    else: norm_score = -(50*(PWM_score-ss_min)/ss_min)
        
    return norm_score

def quant_PWMs(exon_fasta, intron_fasta, PWM_files, soi_file):

	# Parse PWMs
	PWMs = {}
	for f in PWM_files:
		di, ss = f.split('/')[-1].split('.')[0].split('_')
		try: 
			PWMs[di][ss] = { 'pwm': read_PWM_from_file(f), 'min': 0, 'max': 0 }
		except KeyError: 
			PWMs[di] = { ss:  { 'pwm': read_PWM_from_file(f), 'min': 0, 'max': 0 } }

	# Calculate the background scores based on all the 
	# exons and introns supplied
	exons = parse_fasta_file(exon_fasta)
	introns = parse_fasta_file(intron_fasta)

	for t_id in natsorted(exons): 

		if t_id in introns:

			# Get orientation
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

				if dinu in PWMs:

					five_seq = E1[-3:] + I[:10]
					three_seq = I[-14:] + E2[:3]

					try: 
						PWMs[dinu]['5SS'] = update_min_max(five_seq, PWMs[dinu]['5SS'])
						PWMs[dinu]['3SS'] = update_min_max(three_seq, PWMs[dinu]['3SS'])
					except KeyError:
						print "ERROR: Missing either 5' or 3' PWM for {}".format(dinu)

	# Get the scores for the sequences of interest
	# The SOI should have either 5SS or 3SS in 
	# their identifiers.
	soi = {}
	for seq_record in yield_fasta(soi_file):
		seq = seq_record.seq.upper()
		soi[seq_record.id] = { 'seq': seq }
		for ss in ['5SS', '3SS']:
			if ss in seq_record.id.upper():
				soi[seq_record.id]['SS'] = ss
				for dinu in PWMs:
					try:
						score = score_splice_site(seq, PWMs[dinu][ss]['pwm'])
						soi[seq_record.id][dinu] = score
						PWMs[dinu][ss] = update_min_max(seq, PWMs[dinu][ss])
					except KeyError: pass

	print '#ID\tSEQUENCE\tSPLICE_SITE\t{}'.format( '\t'.join([ x for x in natsorted(PWMs) ]) )
	for s in natsorted(soi):
		print '{}\t{}\t{}\t{}'.format(s, soi[s]['seq'], soi[s]['SS'], '\t'.join([ '{}'.format(rescale_score(soi[s][x], PWMs[x][soi[s]['SS']]['min'], PWMs[x][soi[s]['SS']]['max'])) for x in natsorted(PWMs) ]) )

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-e', '--exon-fasta', help="Comma separated list of exon fasta files", required=True)
	parser.add_argument('-i', '--intron-fasta', help="Comma separated list of intron fasta files", required=True)
	parser.add_argument('--PWM', required=True, help="PWM files (output from get_PWMs.py). Files are assumend to look like GT-AG_5SS.txt for example.", nargs='+')
	parser.add_argument('--soi', required=True, help="Fasta files containing the sequences of interest.")
	args = parser.parse_args()

	quant_PWMs(args.exon_fasta, args.intron_fasta, args.PWM, args.soi)