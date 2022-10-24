import sys
import os
#dir_path = os.path.dirname(os.path.dirname(__file__))
#sys.path.insert(0, dir_path)
from pysam import VariantFile
from vcf.functions import getChromLengths, generateWindows
import numpy as np
import argparse

class RunsOfHomozygosity:
	def __init__(self, sample):
		self.sample = sample
		self.chrom = None
		self.start = None
		self.lastSnpPos = None
		self.lastHomPos = None
		self.lastHetPos = None
		self.lastState = None
		self.SnpDensity = None
		self.hetDensity = None
		self.current_length = None
		self.SnpsSinceLastHet = 0
		self.adjacentHets = 0
		self.hom = 0
		self.het = 0
		self.missing = 0
		self.ROH = []
		self.NA = []
	def resetROH(self):
		self.start = None
	def checkSnpDensity(self,pos, min_density=(2,100000)):
		if (self.het + self.hom) / (pos - self.start + 1) > (min_density[0] / min_density[1]):
			return True
		else:
			return False
	def startROH(self, chrom, pos):
		self.start = pos
		self.lastSnpPos = pos
		self.lastHomPos = pos
		self.lastHetPos = None
		self.last_state = "hom"
		self.SnpDensity = 1
		self.hetDensity = 0
		self.SnpsSinceLastHet = 1
		self.adjacentHets = 0
		self.hom = 1
		self.het = 0
		self.missing = 0
	def addHomPos(self,chrom,pos):
		self.lastHomPos = pos
		self.lastSnpPos = pos
		self.lastState = "hom"
		self.adjacentHets = 0
		self.hom += 1
		self.SnpsSinceLastHet += 1
		self.SnpDensity = (self.het + self.hom) / (pos - self.start + 1)
	def addMissing(self,chrom,pos):
		self.missing += 1
	def addMissingTract(self,pos):
		self.NA.append({'chrom':self.chrom,'start':self.lastSnpPos,'end':pos})
		# and reset the ROHs
		self.resetROH()
	def addHetPos(self,chrom,pos):
		self.lastHetPos = pos
		self.lastState = "het"
		self.adjacentHets += 1
		self.SnpsSinceLastHet = 0
		self.SnpDensity = (self.het + self.hom) / (pos - self.start + 1)
	def addROH(self):
		self.ROH.append({
			'chrom': self.chrom,
			'start': self.start,
			'end': self.lastHomPos,
			'length': self.lastHomPos - self.start,
			'density': self.SnpDensity,
			'prophet': self.het / (self.het + self.hom),
			'prophom': self.hom / (self.hom + self.het),
			'propmissing':self.missing / (self.hom + self.het + self.missing)
		})
		self.resetROH()
	
sample = 'PD_0028_Allenopithecus_nigroviridis'
rohs = RunsOfHomozygosity(sample)

for rec in input_vcf.fetch():
	chrom = rec.contig
	pos = rec.pos
	if None in rec.samples[sample]['GT']:
		state = "missing"
	elif set(list(rec.samples[sample]['GT'])) == set([0,1]):
		state = "het"
	elif set(list(rec.samples[sample]['GT'])) in [set([0,0]),set([1,1])]:
		state = "hom"
	# check what we should do with this position
	if rohs.start: #this means that we've got a tract going, and we need to check if we should break it or not
		if state == "hom":
			if pos - rohs.lastSnpPos > max_distance:
				rohs.NA.append({'chrom':chrom, 'start':rohs.lastSnpPos, 'end':pos - 1})
				rohs.startROH(chrom,pos)
			#elif not rohs.checkSnpDensity(pos):
			#	rohs.NA.append({'chrom':chrom,'start':rohs.start,'end':pos - 1})
			#	rohs.startROH(chrom,pos)
			else: # otherwise just add the hompos to the roh
				rohs.addHomPos(chrom,pos)
		elif state == "het": # if it's a het genotype, and we have a tract going, check if there's too many hets in there and break in that case
			if rohs.SnpsSinceLastHet < max_het_per_bp[1]:
				rohs.addROH()
			elif not rohs.checkSnpDensity(pos):
				rohs.NA.append({'chrom':chrom,'start':rohs.start,'end':pos - 1})
				rohs.resetROH()
			else:
				rohs.addHetPos(chrom,pos)
		elif state == "missing":
			rohs.addMissing(chrom,pos)
	else: # if there's no start, change to the correct chromosome and check if we should start a roh
		if rohs.chrom != chrom:
			rohs.chrom = chrom
			rohs.lastSnpPos = 0
		if pos - rohs.lastSnpPos > max_distance:
			rohs.NA.append({'chrom':chrom,'start':rohs.lastSnpPos,'end':pos - 1})
		if state == "hom":
			rohs.startROH(chrom,pos)
		elif state == "het":
			rohs.lastSnpPos = pos
			
max_distance = 1000000
max_het_per_bp = (1,50)

for roh in rohs.ROH:
	if roh['length'] > 100000:
		print(roh)

		if not self.start:

			self.start = rec.pos


#####################################################################################

#set up and parse input arguments

parser = argparse.ArgumentParser(description="Scan inidivudal genomes/genomic regions for runs of homzygosity.")

parser.add_argument('-i', '--input', type=str, help='Input vcf', required=True)

parser.add_argument('-o', '--outprefix', type=str, help='Output prefix.', required=True)

parser.add_argument('--min-density', type=int, help="Minimum SNP density before breaking a ROH, defaults to two per 100 KB (given in count per 100 Kb)", default=2)

parser.add_argument('--max-distance', type=int, help="Maximum distance between adjacent SNPs before breaking ROH, defaults to 1,000,000 bp", default=1000000)

parser.add_argument('--max-adjacent-hetsites', type=int, help="Maximum number of adjacent heterozygous sites to allow before breaking ROH, defaults to 2.", default=2)

parser.add_argument('--max-heterozygosity', type=int, help="Maximum number of heterozygous calls per bp, given as a tuple. defaults to 1 per 50 (1,50)", default = (1,50))

args = parser.parse_args()

# testing
input_vcf = VariantFile("/Users/axeljensen/intro_simulations/bial.pass.AB.gatk.bcftools.NC_041772.1.repfilt.indels_excl.vcf.gz")
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