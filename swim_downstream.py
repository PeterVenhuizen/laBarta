#!/usr/bin/env python

"""
Downstream Salmon analysis.
"""

from __future__ import division
from natsort import natsorted
from fileParser import parse_salmon
import numpy as np
import argparse
import math
import os

from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.sandbox.stats.multicomp import multipletests

def calculate_isoform_ratios(s1, s2, ref_file, labels=['one', 'two'], output_dir="./isoform_ratios/"):

	ref_isoforms = { line.rstrip().split('.')[0]: line.rstrip() for line in open(ref_file) }

	output = []
	for g_id in natsorted(ref_isoforms):
		ref = ref_isoforms[g_id]
		try: 
			for t_id in natsorted(s1[0][g_id].keys()):

				# Add the small value of 0.00001 to avoid zero division errors
				s1_ratios = [ (s1[i][g_id][t_id]+0.00001)/(s1[i][g_id][ref]+0.00001) for i in xrange(0, len(s1)) ]
				s2_ratios = [ (s2[i][g_id][t_id]+0.00001)/(s2[i][g_id][ref]+0.00001) for i in xrange(0, len(s2)) ]
				pval = ttest_ind(s1_ratios, s2_ratios)[1]
				pval = 1 if math.isnan(pval) else pval
				output.append([ g_id, t_id, [ s[g_id][t_id] for s in s1+s2 ], s1_ratios, s2_ratios, pval ])
		except KeyError: pass

	reject, pvals_corrected, alpha_sidak, alpha_bonf = multipletests([ output[i][-1] for i in xrange(len(output)) ], method='fdr_bh')
	with open('{}{}-{}_ratios.txt'.format(output_dir, labels[0], labels[1]), 'w') as fout:
		for i in xrange(len(output)):
			g_id, t_id, TPM, s1_ratios, s2_ratios, pval = output[i]
			#fout.write( '{}\t{}\t{}\t{:.7f}\t{:.7f}\n'.format(t_id, '\t'.join([ '{:.3f}'.format(x) for x in TPM ]), '\t'.join([ '{:.3f}'.format(s) for s in s1_ratios+s2_ratios ]), pval, pvals_corrected[i]) )
			fout.write( '{}\t{}\t{:.7f}\t{:.7f}\n'.format(t_id, '\t'.join([ '{:.3f}'.format(s) for s in s1_ratios+s2_ratios ]), pval, pvals_corrected[i]) )

def find_changes(s1, s2, labels=['one', 'two'], output_dir="./isoform_changes/"):

	output = []
	genes = s1[0].keys()
	for g_id in natsorted(genes):
		t_ids = natsorted(s1[0][g_id].keys())
		for t_id in t_ids:
			#pval = ttest_ind([float(s[g_id][t_id]) for s in s1], [float(s[g_id][t_id]) for s in s2])[1]
			#pval = 1 if math.isnan(pval) else pval
			try: pval = mannwhitneyu([float(s[g_id][t_id]) for s in s1], [float(s[g_id][t_id]) for s in s2])[1]
			except ValueError: pval = 1
			output.append([ t_id, [ s[g_id][t_id] for s in s1+s2 ], pval ])

	reject, pvals_corrected, alpha_sidk, alpha_bonf = multipletests([ output[i][-1] for i in xrange(len(output)) ], method='fdr_bh')
	with open('{}{}-{}_change.txt'.format(output_dir, labels[0], labels[1]), 'w') as fout:
		for i in xrange(len(output)):
			t_id, TPM, pval = output[i]
			fout.write( '{}\t{}\t{:.7f}\t{:.7f}\n'.format(t_id, '\t'.join([ '{:.3f}'.format(x) for x in TPM ]), pval, pvals_corrected[i]) )

def find_isoform_switches(s1, s2, labels=['one', 'two'], output_dir="./isoform_switches/"):

	import subprocess
	import operator

	genes = s1[0].keys()
	for g_id in natsorted(genes):
		t_ids = natsorted(s1[0][g_id].keys())

		# Look for isoform switches
		switches = {}
		N = len(s1)

		for i in xrange(N):

			order1 = sorted(s1[i][g_id].iteritems(), key=operator.itemgetter(1), reverse=True)
			order2 = sorted(s2[i][g_id].iteritems(), key=operator.itemgetter(1), reverse=True)

			names1 = [ o[0] for o in order1 ]
			names2 = [ o[0] for o in order2 ]

			for j, x1 in enumerate(names1):
				x2 = names2[j]
				if x1 != x2:
					one = "before" if names1.index(x1) < names1.index(x2) else "after"
					two = "before" if names2.index(x1) < names2.index(x2) else "after"

					# Only take those switches for which the expression is >= 1 TPM 
					# in at least one condition for each isoform. 
					if one != two and any([ s1[i][g_id][x1] >= 1, s2[i][g_id][x1] >= 1 ]) and any([ s1[i][g_id][x2] >= 1, s2[i][g_id][x2] >= 1 ]):
						try: switches[x1+'-'+x2].append( one+"-"+two )
						except KeyError: switches[x1+'-'+x2] = [ one+'-'+two ]

		# Only output switches which are present in all replicates
		# and make sure that the switches are all in the same direction
		is_switched = False
		printed = set()
		for x in switches:
			if len(switches[x]) == N and len(set(switches[x])) == 1:
				x1, x2, = natsorted(x.split('-'))
				if (x1, x2) not in printed:
					printed.add((x1, x2))
					print '{}\t{}\t{}'.format(g_id, x1, x2)
					is_switched = True

		# Plot the switches (and other isoforms) mean expression
		if len(t_ids) > 1 and is_switched:

			with open('save_dot_plot.R', 'w') as fout:
				fout.write( "#!/usr/bin/env Rscript\n" )
				fout.write( "library(ggplot2)\n" )

				fout.write( "df <- data.frame(\n" )
				fout.write( "\tisoform = rep(c( {} ), each=2),\n".format( ', '.join(['"'+t+'"' for t in t_ids])) )
				fout.write( "\tcondition = rep(c( \"{}\", \"{}\" ), {}),\n".format( labels[0], labels[1], len(t_ids) ) )

				fout.write( "\tTPM_avg = c(\n" )
				fout.write( ',\n'.join([ "\t\t{}, {}".format( np.mean([ s[g_id][t] for s in s1 ]), np.mean([ s[g_id][t] for s in s2 ]) ) for t in t_ids ]) )
				fout.write( "\n\t),\n" )

				fout.write( "\tstdev = c(\n" )
				fout.write( ',\n'.join([ "\t\t{}, {}".format( np.std([ s[g_id][t] for s in s1 ]), np.std([ s[g_id][t] for s in s2 ]) ) for t in t_ids ]) )
				fout.write( "\n\t)\n)")

				fout.write( "\ndf$condition <- factor(df$condition, levels=c(\"{}\", \"{}\"))".format( labels[0], labels[1] ) )
				fout.write( "\nlimits <- aes(ymax = TPM_avg + stdev, ymin = TPM_avg - stdev, colour=isoform)")
				fout.write( "\npng(file=\"{}/{}.png\", width=450, height=600)".format(output_dir, g_id) )
				fout.write( "\nggplot(data=df, aes(x=condition, y=TPM_avg, group=isoform)) +" )
				fout.write( "\n\tggtitle(\"{}\") + theme(plot.title = element_text(hjust = 0.5)) +".format(g_id) )
				fout.write( "\n\tgeom_point(aes(colour=isoform), size=3) +" )
				fout.write( "\n\tgeom_errorbar(limits, width=0.05, size=0.35) +" )
				fout.write( "\n\tgeom_line(size=1, alpha=0.5, aes(colour=isoform)) +" )
				fout.write( "\n\tlabs(x=\"\", y=\"TPM mean\", fill=\"\")")

			subprocess.call("Rscript save_dot_plot.R", shell=True)

	subprocess.call("rm save_dot_plot.R", shell=True)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	subparsers = parser.add_subparsers(dest='command', help='sub-command help')

	parser_a = subparsers.add_parser('get-ratios', help="Calculate the isoform ratios.")
	parser_a.add_argument('-s1', required=True, nargs='+', help="Condition one Salmon file(s).")
	parser_a.add_argument('-s2', required=True, nargs='+', help="Condition two Salmon file(s).")
	parser_a.add_argument('-r', '--ref', required=True, help="List of representative gene model identifiers.")
	parser_a.add_argument('-l', '--labels', nargs='+', default=['one', 'two'], help="Condition labels.")
	parser_a.add_argument('-o', '--output-dir', default="./isoform_ratios/", help="Output directory.")

	parser_b = subparsers.add_parser('get-switches', help="Look for isoform switches.")
	parser_b.add_argument('-s1', required=True, nargs='+', help="Condition one Salmon file(s).")
	parser_b.add_argument('-s2', required=True, nargs='+', help="Condition two Salmon file(s).")
	parser_b.add_argument('-l', '--labels', nargs='+', default=['one', 'two'], help="Condition labels.")
	parser_b.add_argument('-o', '--output-dir', default="./isoform_switches/", help="Output directory.")

	parser_c = subparsers.add_parser('get-changes', help="Look for isoform expression changes.")
	parser_c.add_argument('-s1', required=True, nargs='+', help="Condition one Salmon file(s).")
	parser_c.add_argument('-s2', required=True, nargs='+', help="Condition two Salmon file(s).")
	parser_c.add_argument('-l', '--labels', nargs='+', default=['one', 'two'], help="Condition labels.")
	parser_c.add_argument('-o', '--output-dir', default="./isoform_switches/", help="Output directory.")

	args = parser.parse_args()

	if len(args.command):

		if not os.path.exists(args.output_dir): os.makedirs(args.output_dir)
		s1 = [ parse_salmon(f) for f in args.s1 ]
		s2 = [ parse_salmon(f) for f in args.s2 ]

	if args.command == "get-ratios":
		calculate_isoform_ratios(s1, s2, args.ref, args.labels, args.output_dir+'/')
	elif args.command == "get-switches":
		find_isoform_switches(s1, s2, args.labels, args.output_dir+'/')
	elif args.command == "get-changes":
		find_changes(s1, s2, args.labels, args.output_dir+'/')