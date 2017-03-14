#!/usr/bin/env python

"""
Generate connected dot plots for each gene
and try to identify isoform switches.

To-do: 
	* Output all switched isoforms to stdout
	* Plot the mean values
	* Add error bars
"""

from __future__ import division
from collections import Counter
from natsort import natsorted
from fileParser import parse_salmon
import numpy as np
import subprocess
import argparse
import os
import operator

def all_switched(g_id, s1, s2):

	switches = Counter()
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

				if one != two and all([ s1[i][g_id][x1] >= 1, s1[i][g_id] >= 1, s2[i][g_id] >= 1, s2[i][g_id] >= 1 ]):
					switches[x1+'-'+x2] += 1

	is_switched = False
	printed = set()
	for x in switches:
		if switches[x] == N:
			x1, x2 = natsorted(x.split('-'))
			if (x1, x2) not in printed:
				printed.add((x1, x2))
				print '{}\t{}\t{}'.format(g_id, x1, x2)
				is_switched = True

	return is_switched

def find_isoform_switches(s1_files, s2_files, selection_file, output_dir):
	
	s1 = [ parse_salmon(f) for f in s1_files ]
	s2 = [ parse_salmon(f) for f in s2_files ]

	if not os.path.exists(output_dir): os.makedirs(output_dir)

	genes = s1[0].keys() if selection_file == None else [ line.rstrip() for line in open(selection_file) ]
	for g_id in natsorted(genes):
		t_ids = natsorted(s1[0][g_id].keys())

		#if len(t_ids) > 1 and all([ is_switched(s1[i][g_id], s2[i][g_id]) for i in xrange(0, len(s1)) ]):
		if len(t_ids) > 1 and all_switched(g_id, s1, s2):

			with open('save_dot_plot.R', 'w') as fout:
				fout.write( "#!/usr/bin/env Rscript\n" )
				fout.write( "library(ggplot2)\n" )

				fout.write( "df <- data.frame(\n" )
				fout.write( "\tisoform = rep(c( {} ), each=2),\n".format( ', '.join(['"'+t+'"' for t in t_ids])) )
				fout.write( "\tcondition = rep(c( \"ST\", \"7-1-1\" ), {}),\n".format( len(t_ids) ) )

				fout.write( "\tTPM_avg = c(\n" )
				fout.write( ',\n'.join([ "\t\t{}, {}".format( np.mean([ s[g_id][t] for s in s1 ]), np.mean([ s[g_id][t] for s in s2 ]) ) for t in t_ids ]) )
				fout.write( "\n\t),\n" )

				fout.write( "\tstdev = c(\n" )
				fout.write( ',\n'.join([ "\t\t{}, {}".format( np.std([ s[g_id][t] for s in s1 ]), np.std([ s[g_id][t] for s in s2 ]) ) for t in t_ids ]) )
				fout.write( "\n\t)\n)")

				fout.write( "\ndf$condition <- factor(df$condition, levels=c(\"ST\", \"7-1-1\"))" )
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
	parser.add_argument('-s1', required=True, nargs='+', help="Condition one Salmon files.")
	parser.add_argument('-s2', required=True, nargs='+', help="Condition two Salmon files.")
	parser.add_argument('--labels', nargs='+', default=['one', 'two'], help="Condition labels.")
	parser.add_argument('-o', '--output-dir', default="./", help="Output directory.")
	parser.add_argument('--selection', help="Gene selection file. Gene identifier should be in the first column.")
	args = parser.parse_args()

	find_isoform_switches(args.s1, args.s2, args.selection, args.output_dir)