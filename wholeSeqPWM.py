#!/usr/bin/env python

"""
Whole sequence PWM scoring
"""

from __future__ import division
from natsort import natsorted
from fileParser import yield_fasta
from math import log
import argparse

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
