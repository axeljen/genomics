import sys
import os
dir_path = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, dir_path)
from pysam import VariantFile
from vcf.functions import getChromLengths, generateWindows
import numpy as np
import argparse

def categorizeHetHomMissing(sample, record):
	if None in record.samples[sample]['GT']:
		return "missing"
	else:
		if set(record.samples[sample]['GT']) == set([0,1]):
			return "het"
		elif set(record.samples[sample]['GT']) == set([0,0]) or set(record.samples[sample]['GT']) == set([1,1]):
			return "hom"
		else:
			return "missing"


#####################################################################################

#set up and parse input arguments

parser = argparse.ArgumentParser(description="Calculate per-sample heterozygosity from vcf-data. If --window-size is given, the analyses will be performed in windows along the genome. Otherwise, genomewide heterozygosity is reported.")

parser.add_argument('-i', '--input', type=str, help='Input vcf', required=True)

parser.add_argument('-o', '--outprefix', type=str, help='Output prefix.', required=True)

parser.add_argument('-w', '--window-size', type=int, help='Size of windows to calculate heterozygosity in.', default=None)

parser.add_argument('-s', '--step-size', type=int, help='Step between window starting points.')

parser.add_argument('--samples', type=str, help="File with samples to include, one per line. Defaults to all samples in vcf.")

parser.add_argument('--biallelic-data', action="store_true", help="If specified, heterozygosity will be calculated as heterozygous sites / window size. Otherwise it will be heterozygous sites / total sites.")

parser.add_argument('--min-sites', type=int, help="Minimum number of called sites in a window for heterozygosity to be calculated.", default=1)

parser.add_argument('--output-format', type=str, help="Can be used in conjunction with -w, to specify that the output should be written in wide format rather than the default long (which is easier for plotting e.g. in ggplot", default="long")

args = parser.parse_args()

# testing
#input_vcf = VariantFile("bial.pass.AB.gatk.bcftools.NC_041772.1.repfilt.indels_excl.vcf.gz")
#min_sites = 1
#window_step = 100000
#samples = list(input_vcf.header.samples)
#outprefix = "./test"
# parse the input vcf with pysam
input_vcf = VariantFile(args.input)
# check if samples is in args, otherwise just take all the samples from the vcf
if args.samples:
	if os.isfile(args.samples):
		samples = []
		with open(args.samples) as f:
			for line in f.readlines():
				if not line.rstrip() == "":
					samples.append(line.rstrip()) 
	else:
		samples = args.samples.split(",")
# prep output prefix
outprefix = args.outprefix
# check for window sizes given as inputs 
if args.window_size:
	window_size = args.window_size
	window_analyses = True
else:
	window_size = 100000
	window_analyses = False

if args.step_size:
	window_step = args.step_size
else:
	window_step = window_size

#if no samples are given as input, get list of samples from vcf file
if args.samples:
	samples = args.samples 
else:
	samples = list(input_vcf.header.samples)

#and check whether only biallelic data is given as input
if args.biallelic_data:
	biallelic = True
else:
	biallelic = False

min_sites = args.min_sites
######################################################################################

# If this is a windowed analyses, we'll generate the windows and run this analysis:
if window_size:
	chromlengths = getChromLengths(input_vcf)
	windows = generateWindows(chromlengths, window_size, window_step)

# list for storing window results
results = []

#some checkpoints for printing progress
current_pos = 0
next_checkpoint = current_pos + 10000000
variants = 0
#next iterate over all sites and samples and calculate heterozygosity
with input_vcf as vcf:
	for i, window in enumerate(windows):
		if current_pos == 0:
			print("Calculating heterozygosity for {} samples in windows of {} bp.".format(len(samples),window_size))
		else:
			if current_pos >= next_checkpoint:
				print("Processed " + str(variants) + " variants (" + str(next_checkpoint) + " bp).")
				next_checkpoint = current_pos + 10000000
		window_results = {sample: {'hom':0,'het':0, 'missing':0,} for sample in samples}
		for rec in vcf.fetch(window['chrom'], window['start'], window['end']):
			for sample in samples:
				category = categorizeHetHomMissing(sample, rec)
				window_results[sample][category] +=1
			variants += 1
			current_pos = rec.pos
		#calculate heterozygosity also
		for sample in samples:
			if window_results[sample]['het'] + window_results[sample]['hom'] < args.min_sites:
				window_results[sample]['heterozygosity'] = np.nan
			else:
				if biallelic:
					window_results[sample]['heterozygosity'] = window_results[sample]['het'] / window_size
				else:
					window_results[sample]['heterozygosity'] = window_results[sample]['het'] / (window_results[sample]['het'] + window_results[sample]['hom'])
		results.append(
			{'chrom': window['chrom'], 'start':window['start'] +1, 'end': window['end'], 'results': window_results})

print("Done with calculations, summarizing and writing output.")
# we'll summarize the hetvalues for each sample with numpy
hetlists = {sample: [] for sample in samples}

for result in results:
	for sample in samples:
		hetlists[sample].append(result['results'][sample]['heterozygosity'])
summary_dict = {sample: {} for sample in samples}
for sample in samples:
	a = np.array(hetlists[sample])
	summary_dict[sample]['mean'] = np.nanmean(a)
	summary_dict[sample]['median'] = np.nanmedian(a)
	summary_dict[sample]['stdev'] = np.nanstd(a)
	summary_dict[sample]['var'] = np.nanvar(a)
	summary_dict[sample]['min'] = np.nanmin(a)
	summary_dict[sample]['max'] = np.nanmax(a)
	summary_dict[sample]['5th.q'] = np.nanquantile(a,.05)
	summary_dict[sample]['95th.q'] = np.nanquantile(a,.95)
	summary_dict[sample]['windows'] = len(a)
	summary_dict[sample]['missing'] = np.count_nonzero(np.isnan(a))


# then we start by writing the summarized results
sumout = open(os.path.join(outprefix) + "_hetsum.txt", "w")
sumout.write('\t'.join(['sample','het.mean','het.median','het.stdev','het.var','het.min','het.max','het.5th.q','het.95th.q','het.nwindows','het.nmissing']) + "\n")
for sample in samples:
	mean = str(summary_dict[sample]['mean'])
	median = str(summary_dict[sample]['median'])
	stdev = str(summary_dict[sample]['stdev'])
	var = str(summary_dict[sample]['var'])
	hetmin = str(summary_dict[sample]['min'])
	hetmax = str(summary_dict[sample]['max'])
	het_5th = str(summary_dict[sample]['5th.q'])
	het_95th = str(summary_dict[sample]['95th.q'])
	nwindows = str(summary_dict[sample]['windows'])
	nmissing = str(summary_dict[sample]['missing'])
	sumout.write("\t".join([sample,mean,median,stdev,var,hetmin,hetmax,het_5th,het_95th,nwindows,nmissing]) + "\n")

# if we don't care about the windows, exit here.
if not window_analyses:
	sys.exit(0)

# otherwise we continue with writing the output windows
outwindows = open(os.path.join(outprefix) + "_windowhet.txt", "w")
outwindows.write('\t'.join(['sample','chrom','start','end','heterozygosity','hetsites','homsites','missing']) + "\n")
for window in results:
	chrom = str(window['chrom'])
	start = str(window['start'])
	end = str(window['end'])
	for sample in samples:
		heterozygosity = str(window['results'][sample]['heterozygosity'])
		hetsites = str(window['results'][sample]['het'])
		homsites = str(window['results'][sample]['hom'])
		missing = str(window['results'][sample]['missing'])
		line = '\t'.join([sample,chrom,start,end,heterozygosity,hetsites,homsites,missing]) + "\n"
		outwindows.write(line)