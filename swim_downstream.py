#!/usr/bin/env python

"""
Downstream Salmon/Kallisto basic analysis
"""

from __future__ import print_function, division
from natsort import natsorted
import operator
import argparse
import sys

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

def schooling(s1_files, s2_files, groups=['a', 'b'], mode='TPM', ref_file=None):

	""" 
		Group the Salmon files into one output stream.
		Example output line: G_ID\tT_ID\tA_1\tA_2\tA_3\tB_1\tB_2\t_B3
		The values for A and B can take the shape of
		either the "raw" TPM values, the isoform ratios
		based on the gene reference isoform or the splicing 
		index. 
	"""

	from fileParser import parse_salmon

	s1 = [ parse_salmon(f) for f in s1_files ]
	s2 = [ parse_salmon(f) for f in s2_files ]

	genes = natsorted(s1[0].keys())

	# Output header with group information
	print('# Group A: {}\tN: {}'.format(groups[0], len(s1)))
	print('# Group B: {}\tN: {}'.format(groups[1], len(s2)))

	if mode == 'TPM':
		
		for g_id in genes:
			for t_id in s1[0][g_id].keys():
				print( '{}\t{}\t{}'.format( g_id, t_id, '\t'.join([ str(s[g_id][t_id]) for s in s1+s2 ]) ) )

	else: 

		ref_isoforms = { line.rstrip().split('.')[0]: line.rstrip() for line in open(ref_file) }

		# Tiny factor for avoiding zero division errors
		tiny = 0.000000001

		import operator

		for g_id in genes:

			# Get the reference isoform. If no reference isoform
			# is provided in the file, select the isoform with 
			# the highest expression in the first replicate of 
			# condition one. 
			try: ref = ref_isoforms[g_id]
			except KeyError: ref = max(s1[0][g_id].iteritems(), key=operator.itemgetter(1))[0]

			# Check whether the reference isoform is expressed
			if all([ s[g_id][ref] == 0 for s in s1+s2 ]): eprint('WARNING: Reference isoform {} is not expressed!'.format(ref))

			if mode == 'ratio':
				
				for t_id in s1[0][g_id].keys():
					print('{}\t{}\t{}'.format(g_id, t_id, '\t'.join([ str((s[g_id][t_id]+tiny)/(s[g_id][ref]+tiny)) if s[g_id][ref] else 'NaN' for s in s1+s2 ]) ))

			elif mode == 'splicing_index':
				
				print('{}\t{}\t{}'.format(g_id, ref, '\t'.join([ str((s[g_id][ref]+tiny) / sum([ s[g_id][t_id]+tiny for t_id in s[g_id].keys() ])) if s[g_id][ref] else 'NaN' for s in s1+s2 ]) ))

def get_schooled(school):

	import re
	group_names, group_sizes = [], []

	data = {}
	for line in school:
		if line.startswith('#'):
			name, size = re.search('Group [A|B]: (.*?)\tN: (\d+)', line).groups()
			group_names.append(name)
			group_sizes.append(int(size))
		else:
			cols = line.rstrip().split('\t')
			g_id, t_id = cols[:2]
			try: data[g_id][t_id] = [ float(x) for x in cols[2:] ]
			except KeyError: data[g_id] = { t_id: [ float(x) for x in cols[2:] ] }	

	return data, group_names, group_sizes

def fillet(school, ttest=True, perm_test=True, no_FDR=False):

	"""
		Perform statistical tests on the schooled Salmon data.
	"""

	from stat_tests import mc_perm_test, fdr_correction
	from scipy.stats import ttest_ind

	# Parse school
	data, group_names, group_sizes = get_schooled(school)

	results = {}
	for g_id in data:
		for t_id in data[g_id]:
			if ttest: results[t_id] = { 'ttest': ttest_ind(data[g_id][t_id][:group_sizes[0]], data[g_id][t_id][group_sizes[0]:])[1] }
			if perm_test: 
				try: results[t_id]['perm_test'] = mc_perm_test(data[g_id][t_id][:group_sizes[0]], data[g_id][t_id][group_sizes[0]:])
				except KeyError: results[t_id] = { 'perm_test': mc_perm_test(data[g_id][t_id][:group_sizes[0]], data[g_id][t_id][group_sizes[0]:]) }

	trs = natsorted(results.keys())

	# Do FDR correction for the t-test p-values
	if not no_FDR:
		pval_corrected = fdr_correction([ results[t_id]['ttest'] for t_id in trs ])
		for i, t_id in enumerate(trs):
			results[t_id]['qval'] = pval_corrected[i]

	# Output
	tests = [ 'ttest', 'qval', 'perm_test' ]
	tests = [ tests[i] for i, t in enumerate(tests) if [ttest, not no_FDR, perm_test][i] ]
	header = "#GENE_ID\tTRANSCRIPT_ID\t{}\t{}\t{}".format( 
		'\t'.join([ '{}-{}'.format(group_names[0], i+1) for i in xrange(group_sizes[0]) ]),
		'\t'.join([ '{}-{}'.format(group_names[1], i+1) for i in xrange(group_sizes[1]) ]),
		'\t'.join([ t.upper() for t in tests ])
	)
	print(header)
	for t_id in trs:
		g_id = t_id.split('.')[0]
		print( '{}\t{}\t{}\t{}'.format(g_id, t_id, 
			'\t'.join([ '{:.3f}'.format(x) for x in data[g_id][t_id] ]),
			'\t'.join([ '{:.7f}'.format(results[t_id][x]) for x in tests ])
			) 
		)

def order_isoforms(gene_data, i):

	subset = { t_id: gene_data[t_id][i] for t_id in gene_data }
	order = sorted(subset.iteritems(), key=operator.itemgetter(1), reverse=True)
	return [ o[0] for o in order ]

def find_isoform_switches(school, output_dir='./isoform_switches/'):

	import numpy as np
	import os
	import subprocess

	if not os.path.exists(output_dir): os.makedirs(output_dir)

	# Parse school
	data, labels, group_sizes = get_schooled(school)

	for g_id in natsorted(data):
		t_ids = data[g_id].keys()

		# Look for isoform switches
		if len(t_ids) > 1:
			
			switches = {}
			N = group_sizes[0]

			for i in xrange(N):

				order1 = order_isoforms(data[g_id], i)
				order2 = order_isoforms(data[g_id], i+N)

				for j, x1 in enumerate(order1):
					x2 = order2[j]
					if x1 != x2:
						one = "before" if order1.index(x1) < order1.index(x2) else "after"
						two = "before" if order2.index(x1) < order2.index(x2) else "after"

						# Only take those switches for which the expression is >= 1 TPM
						# in at least one condition for each isoform
						if one != two and any([ data[g_id][x1][i] >= 1, data[g_id][x1][i+N] >= 1 ]) and any([ data[g_id][x2][i] >= 1, data[g_id][x2][i+N] >= 1 ]):
							try: switches[x1+'-'+x2].append( one+'-'+two )
							except KeyError: switches[x1+'-'+x2] = [ one+'-'+two ]

			# Only output switches which are present and going in the 
			# same direction in all replicates. 
			is_switched, printed = False, set()
			for x in switches:
				if len(switches[x]) == N and len(set(switches[x])) == 1:
					x1, x2 = natsorted(x.split('-'))
					if (x1, x2) not in printed:
						printed.add((x1, x2))
						print('{}\t{}\t{}'.format(g_id, x1, x2))
						is_switched = True

			if is_switched:

				with open('save_dot_plot.R', 'w') as fout:
					fout.write( "#!/usr/bin/env Rscript\n" )
					fout.write( "library(ggplot2)\n" )

					fout.write( "df <- data.frame(\n" )
					fout.write( "\tisoform = rep(c( {} ), each=2),\n".format( ', '.join([ '"'+t+'"' for t in t_ids ]) ) )
					fout.write( "\tcondition = rep(c( \"{}\", \"{}\" ), {}),\n".format( labels[0], labels[1], len(t_ids) ) )

					fout.write( "\tTPM_avg = c(\n" )
					fout.write( ',\n'.join([ "\t\t{}, {}".format( np.mean(data[g_id][t][:N]), np.mean(data[g_id][t][N:]) ) for t in t_ids ]) )
					fout.write( "\n\t),\n" )

					fout.write( "\tstdev = c(\n" )
					fout.write( ',\n'.join([ "\t\t{}, {}".format( np.std(data[g_id][t][:N]), np.std(data[g_id][t][N:]) ) for t in t_ids ]) )
					fout.write( "\n\t)\n)" )

					fout.write( "\ndf$condition <- factor(df$condition, levels=c(\"{}\", \"{}\"))".format( labels[0], labels[1] ) )
					fout.write( "\nlimits <- aes(ymax = TPM_avg + stdev, ymin = TPM_avg - stdev, colour = isoform)" )
					fout.write( "\npng(file=\"{}/{}.png\", width=450, height=600)".format(output_dir, g_id) )
					fout.write( "\nggplot(data=df, aes(x=condition, y=TPM_avg, group=isoform)) +" )
					fout.write( "\n\tggtitle(\"{}\") + theme(plot.title = element_text(hjust = 0.5)) +".format(g_id) )
					fout.write( "\n\tgeom_point(aes(colour=isoform), size=3) +" )
					fout.write( "\n\tgeom_errorbar(limits, width=0.05, size=0.35) +" )
					fout.write( "\n\tgeom_line(size=1, alpha=0.5, aes(colour=isoform)) +" )
					fout.write( "\n\tlabs(x=\"\", y=\"TPM mean\", fill=\"\")" )

				subprocess.call("Rscript save_dot_plot.R", shell=True)
				subprocess.call("rm save_dot_plot.R", shell=True)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	subparsers = parser.add_subparsers(dest='cmd', help="Sub-command help")

	parser_a = subparsers.add_parser('schooling', help="Group Salmon (school of fish) files and prepare for statistical analysis.")
	parser_a.add_argument('-s1', required=True, nargs='+', help="Condition A Salmon file(s).")
	parser_a.add_argument('-s2', required=True, nargs='+', help="Condition B Salmon file(s).")
	parser_a.add_argument('--mode', choices=['TPM', 'ratio', 'splicing_index'], default='TPM', help="Mode of preparing the data. 'TPM' => Output TPM value per transcript, 'ratio' => Output isoform ratio per transcript, and 'splicing_index' => Calculate the splicing index per gene.")
	parser_a.add_argument('-g', '--groups', nargs='+', default=['A', 'B'], help="Group labels/names.")
	parser_a.add_argument('-r', '--ref', help="List of representative gene model identifiers.")

	parser_b = subparsers.add_parser('fillet', help="Perform statistical tests on the school of Salmon.")
	parser_b.add_argument('-s', '--school', type=argparse.FileType('r'), default='-', help="Schooled Salmon file. Output from running swim_downstream.py schooling.")
	parser_b.add_argument('--ttest', action="store_true", help="Perform a t-test on the data. By default the p-values are FDR corrected. This can be disabled by adding --no-FDR.")
	parser_b.add_argument('--no-FDR', action="store_true", help="Don't do FDR correction for the t-test.")
	parser_b.add_argument('--perm-test', action="store_true", help="Perform a permutation test on the data.")

	parser_c = subparsers.add_parser('switch', help="Look for isoform switches in the school of Salmon.")
	parser_c.add_argument('-s', '--school', type=argparse.FileType('r'), default='-', help="Schooled Salmon file. Output from running swim_downstream.py schooling.")
	parser_c.add_argument('-o', '--output-dir', default='./isoform_switches', help="Output directory.")

	args = parser.parse_args()

	if args.cmd == "schooling":

		# Check whether a reference file was supplied when the
		# ratio or splicing_index option are selected.
		if args.mode in ['ratio', 'splicing_index'] and args.ref == None:
			parser.error("--mode {ratio|splicing_index} requires --ref.")

		schooling(args.s1, args.s2, args.groups, args.mode, args.ref)

	elif args.cmd == "fillet":

		# Check whether a stastical test was given
		if all([ not args.ttest, not args.perm_test ]): 
			parser.error("fillet requires either --test and/or --perm-test.")

		fillet(args.school, args.ttest, args.perm_test, args.no_FDR)

	elif args.cmd == "switch":

		find_isoform_switches(args.school, args.output_dir)