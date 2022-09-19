import argparse
import sys
import seqtools.functions as seqtools
from seqtools.functions import FastaFetch, ReadPhylip
import pysam
import vcf.functions as v
import os

def parseInput(inputfile):
	if inputfile.endswith((".fa",".fasta",".fa.gz",".fasta.gz")):
		return FastaFetch(inputfile)
	elif inputfile.endswith((".phy",".phylip")):
		return ReadPhylip(inputfile)
	else:
		print("Unknown input format.")
def subsetSeqDict(seqDict, samples):
	seqDict_f = {}
	for k,v in seqDict.items():
		if k in samples:
			seqDict_f[k] = v
	return seqDict_f

def writeOutput(outfile,seqDict):
	if outfile.endswith((".fa",".fasta",".fa.gz",".fasta.gz")):
		return seqtools.WriteSeqDictToFasta(seqDict, outfile)
	elif outfile.endswith((".phy",".phylip")):
		return seqtools.WriteSeqDictToPhylip(seqDict,outfile)
	else:
		print("Unknown output format.")


parser = argparse.ArgumentParser(description="This tool will convert the records in a VCF file to a multi sequence alignment.")

# input and output
parser.add_argument('-i', '--input', type=str, help='Input alignment', required=True)
parser.add_argument('-o', '--output', type=str, help='Output file in fasta (ending with .fa or .fa.gz) or in phylip (ending with .phy or phy.gz).', required=True)

# samples
parser.add_argument('--samples', type=str, help="File with samples to include, one per line.")

args = parser.parse_args()

# read in the file
seqDict = parseInput(args.input)

# if no samples are given, we'll keep all samples
if not args.samples:
	samples = list(seqDict.keys())
else:
	samples = seqtools.parseSampleFile(args.samples)

# filter based on samples
sf = subsetSeqDict(seqDict, samples)

# write outfile
writeOutput(args.output, sf)