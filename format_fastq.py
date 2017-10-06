#! /usr/bin/env python

import argparse
import gzip
import os

from warnings import warn
from Bio.SeqIO.QualityIO import FastqGeneralIterator

desc = """Should be run from inside a Project_XX directory. Takes a directory of
fastq or fastq.gz files, or Sample_XX directories (or any combination thereof)
and converts any dual-indexed Casava-produced fastq(.gz) or single- or
dual-indexed, Basespace-produced fastq(.gz) files into a format equivalent to
that of the single-indexed, Casava-produced fastq(.gz) files. Gzips files that
are already .gz. 

This script must be run on a directory that only includes files that would go in
a single FOFN (i.e. the files must be divided out by second index and put in
different Project directories) or there is the risk that files with identical
first indices will overwrite each other.

The filename formats are (where X are numbers and N are {A,T,C,G}):
Basespace:				samp_SX_LXXX_RX_XXX.fastq(.gz) 
dual-index Casava:		samp_NNNNNN-NNNNNN_LXXX_RX_XXX.fastq(.gz)
single-index Casava:	samp_NNNNNN_LXXX_RX_XXX.fastq(.gz)

The read ID line formats are
dual-index Basespace:	@M00719:295:000000000-AN59N:1:1101:16589:2209 1:N:0:ACAGTG+ATTGGC
single-index Basespace:	@M00719:295:000000000-AN59N:1:1101:16589:2209 1:N:0:ACAGTG
dual-index Casava:		@MISEQ:295:000000000-AN59N:1:1101:16589:2209 1:N:0:ACAGTGATTGGC
single-index Casava:	@MISEQ:295:000000000-AN59N:1:1101:16589:2209 1:N:0:ACAGTG
"""


parser = argparse.ArgumentParser(description = desc)

hstr = 'By default, the original files are deleted to save disk space. Use '+\
		'this flag to keep the original files.'
parser.add_argument('-k',help = hstr, action='store_true')

hstr = 'The directory containing the files you wish to fix. Defaults to the '+\
		'working directory. This directory is expected to contain either'+\
		'.fastq(.gz) files or Sample_XX directories and the script will '+\
		'ignore	all other files.'
parser.add_argument('-i',help = hstr, type=str, default='.')

hstr = 'Verbose. Without this flag, only warnings will be shown.'
parser.add_argument('-v',help=hstr,action='store_true')

def basespace_check(fname):
	'''A function that takes a file name string and checks if the file was
	produced by basespace or casava. Returns True or False, respectively.'''

	nlst = fname.split('_')
	samp = nlst[-4]

	return (samp.startswith('S'))

def barcode_check_get(string,bspc):
	'''A function that takes either a file name string (if the file was
	basespace) or a read ID line string (if the file was casava), and a switch
	(bspc) indicating which it is. Checks if the file is single- or dual-indexed
	and returns a switch (T if dual) and the first barcode.'''
	
	if bspc:
		# Split on ':' and get the first index a second one exists
		longbar = string.split(':')[-1]
		dual = ('+' in longbar)
		first = longbar.split('+')[0]

	else:
		# Split on '_' and get the first index if a second one exists
		longbar = string.split('_')[-4]
		dual =  ('-' in longbar)
		first = longbar.split('-')[0]

	return(dual,first)


def update_fname(fname,barcode):
	'''A function that takes a string formatted like a bcl2fastq2 Illumina file
	name (with or without .gz) and and the new barcode to be associated with
	that file and converts it to a single-barcoded casava Illumina fastq file
	name. Returns a string of the new file name.'''

	nlst = fname.split('_')
	nlst[-4] = barcode

	return '_'.join(nlst)

def update_id(recid,first):
	'''A function that takes a recid string (without the @, as produced by
	FastqGeneralIterator) and the first index. Returns the correctly-formatted
	recid string (still without the @).'''

	reclst = recid.split(':')
	reclst[0] = 'MISEQ'		# This does no harm if the file is already casava
	reclst[-1] = first		# Update the barcode

	recid = ':'.join(reclst)
	return(recid)

def open_fastq(fastq):
	'''A function that takes a fastq(.gz) file name, checks if it's gzipped, and
	returns the a switch and the file handle.'''
	
	gz = fastq.endswith('gz')
	if gz:
		handle = gzip.open(fastq)
	else:
		handle = open(fastq)

	return(gz,handle)

def run_fastqs(path,v):
	'''A function that operates inside the 'path' directory. Runs the actual
	second-barcode removal (if necessary), updates the file name, and re-writes
	the files.'''

	# Iterate through all the files in the directory
	allfiles = os.listdir(path)

	for fastq in allfiles:
		
		# Look for Sample_XX directories
		if fastq.startswith('Sample'):
			# Make the new path string and recurse
			npath = path+fastq+'/'
			run_fastqs(npath,v)

		# Skip everything that isn't a fastq or fastq.gz file
		if (not fastq.endswith('.fastq')) and (not fastq.endswith('.fastq.gz')):
			continue

		# Prepend the directory name so it doesn't matter where we are
		fastq = path+fastq
		
		if v:
			print('PROCESSING: %s' % fastq)

		# Check for gzip and create a file handle
		gz,handle = open_fastq(fastq)

		# Check if the file is basespace or casava
		bspc = basespace_check(fastq)

		# If the file is casava, check its indexing
		if not bspc:
			dual,first = barcode_check_get(fastq,bspc)

			if (not dual) and v:	# this is a single-indexed Casava file. Move on.
				# Skip
				msg = '%s is already single-indexed Casava. Nothing to do.' %\
						fastq
				print(msg)
				continue

		else:
			# If it's basespace, we'll need the first line of the file to get
			# the barcode. Then reset the pointer to the top of the file.
			fline = handle.next().strip()
			gz,handle = open_fastq(fastq)
			
			dual,first = barcode_check_get(fline,bspc)

		# Create a list to store the sequence strings
		seqs = []

		# Iterate with the FastqGeneralIterator
		recs = FastqGeneralIterator(handle)
		for recid,seq,qual in recs:
	
			# Construct the new record string
			recid = update_id(recid,first)
			rec = '@%s\n%s\n+\n%s\n' % (recid,seq,qual)
	
			# Store the records
			seqs.append(rec)

		# Create the new filename with the barcode
		newfq = update_fname(fastq,first)

		# Check if the file exists
		if newfq.split('/')[-1] in allfiles:
			msg = 'WARNING: %s exists and is being overwritten' % newfq
			print(msg)

		# Open the output file for writing
		if gz:
			outf = gzip.open(newfq,'w')
		else:
			outf = open(newfq,'w')
	
		# Write the new file
		for rec in seqs:
			outf.write(rec)

		handle.close()
		outf.close()

		# Delete the old file
		if not keep:
			if fastq != newfq:
				os.remove(fastq)

	return

if __name__ == '__main__':

	args = parser.parse_args()
	keep = args.k
	here = args.i
	if not here.endswith('/'):
		here = here+'/'
	v = args.v

	# Which level are we operating at?
	run_fastqs(here,v)
