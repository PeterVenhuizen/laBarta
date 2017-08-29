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

def get_intron_info(exon_fasta, intron_fasta, PWM_files, exitron_file="", SUPPA_files=[]):

	# Parse PWMs
	PWMs = {}
	for f in PWM_files:
		file_info = f.split('/')[-1].split('.')[0].split('_')
		di, ss = file_info[:2]
		try: 
			PWMs[ss]['_'.join(file_info)] = { 'pwm': read_PWM_from_file(f), 'min': 0, 'max': 0 }
		except KeyError:
			PWMs[ss] = { '_'.join(file_info): { 'pwm': read_PWM_from_file(f), 'min': 0, 'max': 0 } }

	five_keys = natsorted(PWMs['5SS'].keys())
	three_keys = natsorted(PWMs['3SS'].keys())

	# Parse exons and introns
	exons = parse_fasta_file(exon_fasta)
	introns = parse_fasta_file(intron_fasta)

	# Parse SUPPA events if supplied
	suppa = {}

	# Parse exitrons
	for line in open(exitron_file):
		c, s, e = line.split('\t')[:3]
		suppa['{}:{}-{}'.format(c, s, int(e)-1)] = ["EI"]

	for f in SUPPA_files:
		for line in open(f):
			seqname, gene_id, event_id, inclusion_transcripts, total_transcripts = line.rstrip().split('\t')

			events = []
		
			try:
				event_type = event_id.split(';')[1].split(':')[0]

				if event_type == 'A3':
					event_type, seqname, e1_s2, e1_s3, strand = event_id.split(';')[1].split(':')

					# Short
					e1, s2 = e1_s2.split('-')
					events.append(['{}:{}-{}'.format(seqname, e1, int(s2)-1), 'A3'])

					# Long
					e1, s3 = e1_s3.split('-')
					events.append(['{}:{}-{}'.format(seqname, e1, int(s3)-1), 'A3'])

				elif event_type == 'A5':
					event_type, seqname, e2_s3, e1_s3, strand = event_id.split(';')[1].split(':')

					# Short
					e2, s3 = e2_s3.split('-')
					events.append(['{}:{}-{}'.format(seqname, e2, int(s3)-1), 'A5'])

					# Long
					e1, s3 = e1_s3.split('-')
					events.append(['{}:{}-{}'.format(seqname, e1, int(s3)-1), 'A5'])

				elif event_type == 'RI':
					event_type, seqname, s1, e1_s2, e2, strand = event_id.split(';')[1].split(':')
					e1, s2 = e1_s2.split('-')
					coord = '{}:{}-{}'.format(seqname, e1, int(s2)-1)
					try: 
						if 'EI' not in suppa[coord]:
							events.append([coord, 'IR'])
					except KeyError: 
						events.append([coord, 'IR'])

				elif event_type == 'SE':
					event_type, seqname, e1_s2, e2_s3, strand = event_id.split(';')[1].split(':')
					e1, s2 = e1_s2.split('-')
					e2, s3 = e2_s3.split('-')
					events.append(['{}:{}-{}'.format(seqname, e1, int(s3)-1), 'ES'])

				elif event_type == 'MX':
					event_type, seqname, e1_s2, e2_s4, e1_s3, e3_s4, strand = event_id.split(';')[1].split(':')
					e1, s3 = e1_s3.split('-')
					events.append(['{}:{}-{}'.format(seqname, e1, int(s3)-1), 'MX'])

					e2, s4 = e2_s4.split('-')
					events.append(['{}:{}-{}'.format(seqname, e2, int(s4)-1), 'MX'])

			except IndexError: pass

			for e in events:
				try: suppa[e[0]].append(e[1])
				except KeyError: suppa[e[0]] = [e[1]]

	# Collect intron information
	info = {}
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

				# Get the splice site sequences
				five_seq = E1[-3:] + I[:10]
				three_seq = I[-14:] + E2[:3]
				dinu = I[:2]+'-'+I[-2:]

				if t_id not in info: info[t_id] = []
				info[t_id].append({
					'chr': c,
					'coord': x,
					'strand': '+' if isFwd else '-',
					'dinu': dinu,
					'length': len(I),
					'GC': (I.count('G')+I.count('C'))/len(I),
					'5SS': five_seq,
					'5SS_scores': [],
					'3SS': three_seq,
					'3SS_scores': [],
					'suppa': ','.join(natsorted(set(suppa[c+':'+x]))) if c+':'+x in suppa else "."
				})

				for k in five_keys:
					info[t_id][-1]['5SS_scores'].append(score_splice_site(five_seq, PWMs['5SS'][k]['pwm']))
					PWMs['5SS'][k] = update_min_max(five_seq, PWMs['5SS'][k])

				for k in three_keys:
					info[t_id][-1]['3SS_scores'].append(score_splice_site(three_seq, PWMs['3SS'][k]['pwm']))
					PWMs['3SS'][k] = update_min_max(three_seq, PWMs['3SS'][k])

	print '#TRANSCRIPT_ID\tCHR\tCOORDINATES\tSTRAND\tDINUCLEOTIDES\tAS\tLENGTH\tGC_CONTENT\t5SS_SEQ\t5SS_INTRON_TYPE\t5SS_BEST_SCORE\t3SS_SEQ\t3SS_INTRON_TYPE\t3SS_BEST_SCORE\tINTRON_LABEL'
	#print '#TRANSCRIPT_ID\tCHR\tCOORDINATES\tDINUCLEOTIDES\t5SS_INTRON_TYPE\t'+'\t'.join(five_keys)+'\t3SS_INTRON_TYPE\t'+'\t'.join(three_keys)+'\tINTRON_LABEL'
	for t_id in natsorted(info):
		for i in xrange(len(info[t_id])):

			# Label introns
			scaled_5 = [ rescale_score(info[t_id][i]['5SS_scores'][j], PWMs['5SS'][k]['min'], PWMs['5SS'][k]['max']) for j, k in enumerate(five_keys) ]
			label_5 = five_keys[scaled_5.index(max(scaled_5))] if max(scaled_5) >= 65 else 'LOWSCORE_5SS'

			scaled_3 = [ rescale_score(info[t_id][i]['3SS_scores'][j], PWMs['3SS'][k]['min'], PWMs['3SS'][k]['max']) for j, k in enumerate(three_keys) ]
			label_3 = three_keys[scaled_3.index(max(scaled_3))] if max(scaled_3) >= 65 else 'LOWSCORE_3SS'

			max5 = five_keys[scaled_5.index(max(scaled_5))]
			max3 = three_keys[scaled_3.index(max(scaled_3))]

			if max(scaled_5) >= 65:

				if 'U12' in max5 and max(scaled_5) >= 75 and 'U12' in max3 and max(scaled_3) >= 65:
					intron_label = max5.replace('_5SS', '')
				elif 'U12' in max5 and max(scaled_5) < 75:

					if 'GT-AG' in max5: 
						intron_label = max5.replace('_5SS', '')

					# Get 2nd largest and check if it still is over 65
					elif sorted(scaled_5)[-2] >= 65:
						intron_label = five_keys[scaled_5.index(sorted(scaled_5)[-2])].replace('_5SS', '')
				else: 
					intron_label = max5.replace('_5SS', '')

			elif max(scaled_3) >= 65:
				
				dinu_U2 = '{}_3SS_U2'.format(info[t_id][i]['dinu'])
				try:
					if scaled_3[three_keys.index(dinu_U2)] >= 65: 
						intron_label = dinu_U2.replace('_3SS', '')
					else: 
						intron_label = "LOWSCORE"
				except ValueError: 
					intron_label = "LOWSCORE"

			else:
				intron_label = "LOWSCORE"

			"""
			print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
				t_id, 
				info[t_id][i]['chr'],
				info[t_id][i]['coord'],
				info[t_id][i]['dinu'],
				label_5,
				'\t'.join([ '{:.3f}'.format(x) for x in scaled_5 ]),
				label_3,
				'\t'.join([ '{:.3f}'.format(x) for x in scaled_3 ]),
				intron_label
			)
			"""

			print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\t{:.3f}\t{}\t{}\t{:.3f}\t{}'.format(
				t_id, 
				info[t_id][i]['chr'],
				info[t_id][i]['coord'],
				info[t_id][i]['strand'],
				info[t_id][i]['dinu'],
				info[t_id][i]['suppa'],
				info[t_id][i]['length'],
				info[t_id][i]['GC'],
				info[t_id][i]['5SS'],
				label_5,
				max(scaled_5),
				info[t_id][i]['3SS'],
				label_3,
				max(scaled_3),
				intron_label
			)


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-e', '--exon-fasta', help="Comma separated list of exon fasta files", required=True)
	parser.add_argument('-i', '--intron-fasta', help="Comma separated list of intron fasta files", required=True)
	parser.add_argument('--PWM', required=True, help="PWM files (output from get_PWMs.py). Files are assumed to look like GT-AG_5SS.txt for example.", nargs='+')
	parser.add_argument('--exitrons', required=True, help="Exitrons file.")
	parser.add_argument('--SUPPA', nargs='+', default=[], help="SUPPA ioe strict files.")
	args = parser.parse_args()

	get_intron_info(args.exon_fasta, args.intron_fasta, args.PWM, args.exitrons, args.SUPPA)