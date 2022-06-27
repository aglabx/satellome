#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 29.01.2013
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

import sys
from trseeker.tools.trf_tools import trf_search_in_dir_parallel

if __name__ == '__main__':
	if len(sys.argv) != 5:
		print "Usage: name.py input_folder output_folder mask threads"
		exit(0)
	else:
		input_folder = sys.argv[1]
		output_folder = sys.argv[2]
		mask = sys.argv[3]
		threads = sys.argv[4]
	threads = 40
	trf_search_in_dir_parallel(input_folder, verbose=True, file_suffix=mask, output_folder=output_folder, threads=threads)