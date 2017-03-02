#!/usr/bin/env python

"""
Generate connected dot plots for each gene
and try to identify isoform switches.
"""

from __future__ import division
from natsort import natsorted
from fileParser import parse_salmon
import subprocess
import argparse
import os
import operator

def is_switched(trs1, trs2):
	
	order1 = sorted(trs1.iteritems(), key=operator.itemgetter(1), reverse=True)
	order2 = sorted(trs2.iteritems(), key=operator.itemgetter(1), reverse=True)

	names1 = [ o[0] for o in order1 ]
	names2 = [ o[0] for o in order2 ]

	for i, x1 in enumerate(names1):
		x2 = names2[i]
		if x1 != x2:
			one = "before" if names1.index(x1) < names1.index(x2) else "after"
			two = "before" if names2.index(x1) < names2.index(x2) else "after"

			if one != two:
				if all([ trs1[x1] >= 1, trs1[x2] >= 1, trs2[x1] >= 1, trs2[x2] >= 1 ]):
					return True

	return False

def plotify(s1_files, s2_files, selection_file, output_dir):
	
	s1 = [ parse_salmon(f) for f in s1_files ]
	s2 = [ parse_salmon(f) for f in s2_files ]

	if not os.path.exists(output_dir): os.makedirs(output_dir)

	genes = s1[0].keys() if selection_file == None else [ line.rstrip() for line in open(selection_file) ]
	for g_id in natsorted(genes):
		t_ids = natsorted(s1[0][g_id].keys())

		if len(t_ids) > 1 and all([ is_switched(s1[i][g_id], s2[i][g_id]) for i in xrange(0, len(s1)) ]):
			print g_id

			with open('save_dot_plot.R', 'w') as fout:
				fout.write( "#!/usr/bin/env Rscript\n" )
				fout.write( "library(ggplot2)\n" )

				fout.write( "df <- data.frame(\n" )
				fout.write( "\tisoform = rep(c( {} ), each={}),\n".format( ', '.join(['"'+t+'"' for t in t_ids]), len(s1)+len(s2)) )
				fout.write( "\tcondition = rep(c( \"ST\", \"7-1-1\" ), each={}),\n".format( len(s1) ) )
				
				fout.write( "\tfullname = c(\n" )
				fout.write(  ',\n'.join([ "\t\trep(c({}), 2)".format( ', '.join([ '"'+t+"_"+str(r+1)+'"' for r in xrange(len(s1)) ]) ) for t in t_ids ]) )
				fout.write( "\n\t),\n" )

				fout.write( "\tPWM = c(\n" )
				fout.write( ',\n'.join([ "\t\t{}".format( ', '.join([ str(s[g_id][t]) for s in s1+s2 ]) ) for t in t_ids ]))
				fout.write( "\n\t)\n)" )

				fout.write( "\ndf$condition <- factor(df$condition, levels=c(\"ST\", \"7-1-1\"))" )
				fout.write( "\npng(file=\"{}/{}.png\", width=450, height=600)".format(output_dir, g_id) )
				fout.write( "\nggplot(data=df, aes(x=condition, y=PWM, group=fullname)) +" )
				fout.write( "\n\tggtitle(\"{}\") + theme(plot.title = element_text(hjust = 0.5)) +".format(g_id) )
				fout.write( "\n\tgeom_point(aes(colour=isoform), size=3) +" )
				fout.write( "\n\tgeom_line(size=1, alpha=0.5, aes(colour=isoform)) +" )
				fout.write( "\n\tlabs(x=\"\", y=\"PWM\", fill=\"\")")

			subprocess.call("Rscript save_dot_plot.R", shell=True)


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-s1', required=True, nargs='+', help="Condition one Salmon files.")
	parser.add_argument('-s2', required=True, nargs='+', help="Condition two Salmon files.")
	parser.add_argument('--labels', nargs='+', default=['one', 'two'], help="Condition labels.")
	parser.add_argument('-o', '--output-dir', default="./", help="Output directory.")
	parser.add_argument('--selection', help="Gene selection file. Gene identifier should be in the first column.")
	args = parser.parse_args()

	plotify(args.s1, args.s2, args.selection, args.output_dir)