
def parsePopfile(popfile):
	clades = {}
	with open(popfile) as f:
		for line in f.readlines():
			sample = line.split("\t")[0]
			try:
				clade = line.split("\t")[1].rstrip()
				if not clade in list(clades.keys()) and not clade == "":
					clades[clade] = [sample]
				else:
					clades[clade].append(sample)
			except:
				continue
	return clades

def DifferentiallyFixed(pops, seqDict):
	diff = False
	popalleles = {}
	# count missing
	missing = {pop: 0 for pop in pops.keys()}
	for pop,samples in pops.items():
		popalleles[pop] = []
		for sample in samples:
			if not seqDict[sample][0] in popalleles[pop] and seqDict[sample][0] in ['A','T','C','G','a','t','c','g']:
				popalleles[pop].append(seqDict[sample][0])
			elif seqDict[sample][0] in ['M','m','R','r','W','w','S','s','Y','y','K','k']:
				alleles = seqtools.IUPACReverse(seqDict[sample][0])
				for a in alleles:
					if not a in popalleles[pop] and a in ['A','T','C','G','a','t','c','g']:
						popalleles[pop].append(a)
			elif seqDict[sample][0] in ['N','n','-','.']:
				missing[pop] += 1
	# change missing count to fraction
	for pop,i in missing.items():
		missing[pop] = i / len(pops[pop])
	# check if the site is differentially fixed between at least two pops
	sitealleles = []
	for pop,alleles in popalleles.items():
		if len(alleles) == 1:
			if not alleles[0] in sitealleles:
				sitealleles.append(alleles[0])
	if len(set(sitealleles)) > 1:
		diff = True
		for pop in pops.keys():
			if len(popalleles[pop]) == 0:
				popalleles[pop] = 'na'
			elif len(popalleles[pop]) == 1:
				continue
			elif len(popalleles[pop]) > 1:
				popalleles[pop] = list(set(popalleles[pop]))
		return popalleles, missing
	else:
		return False, False

	
import seqtools.functions as seqtools
import argparse, re, pysam
from vcf.functions import get_bases, Haploidize

parser = argparse.ArgumentParser()

#arguments for in and output:
parser.add_argument("-i", "--input", help="Input file, either as alignment in fasta format or as a vcf file.", required=True)
parser.add_argument("-r", "--region", help="Region to analyze as chrom:start-end, compatible only with vcf input.")
parser.add_argument("-o", "--output-prefix", help="Prefix to output file(s).", required=True)

parser.add_argument("-p", "--popfile", help="File with two tab separated columns, sample in the first and clade/population in second. It's ok to have samples without assignment.", required=True)
parser.add_argument("--max-missing", help="Fraction of maximum of missing samples per population to allow for outputting a fixed difference.", default = 1)

parser.add_argument("--verbose-output", help="Added this option to give a more informative output field", action="store_true")
parser.add_argument("--output-vcf", help="With this flag, we'll output the vcf record for sites identified as fixed, for simple parsing to VEP.", action="store_true")

args = parser.parse_args()


#fasta = "/Users/axeljensen/Dropbox/FRANKIES_PROJECT/TFB2M.fa"

# read popfile and parse it
popfile=args.popfile
pops = parsePopfile(popfile)

# loop through each position to look for fixed differences
fixed_diffs = {}
diffs_found = 0

# try to parse the input file format, if failing we will exit 
input_read = False

# max missing from args input
maxmissing = float(args.max_missing)

# if a fasta input is given, take this path
if re.match("^.*(.fa|.fa.gz|.fasta|.fasta.gz)$", args.input):
	input_read = True
	seqDict = seqtools.FastaFetch(args.input_fasta)
	# check length of alignment
	seqlen = len(seqDict[list(seqDict.keys())[0]])
	for base in range(0,seqlen):
		posDict = {sample: seqDict[sample][base] for sample in seqDict.keys()}
		diff, missing = DifferentiallyFixed(pops, posDict)
		if diff:
			for m in missing.values():
				if not m > maxmissing:
					fixed_diffs[base + 1] = diff
					diffs_found +=1

# otherwise, check that it's a vcf and continue this path instead
if re.match("^.*(.vcf|.vcf.gz)$", args.input):
	input_read = True
	samples = []
	for pop,sample_list in pops.items():
		samples = samples + sample_list
	vcf = pysam.VariantFile(args.input)
	if args.region:
		chrom,interval = args.region.split(":")[0],args.region.split(":")[1]
		start,end = int(interval.split("-")[0]),int(interval.split("-")[1])
	else:
		chrom,start,end = None,None,None
	for rec in vcf.fetch(chrom,start,end):
		chrom = rec.chrom
		base = rec.pos
		ref = rec.ref
		alts = ','.join(rec.alts)
		sample_alleles = get_bases(rec,samples)
		posDict = Haploidize(sample_alleles, use_ambiguities=True)
		diff, missing = DifferentiallyFixed(pops,posDict)
		if diff:
			maxmissing_filtered = False
			for m in missing.values():
				if m > maxmissing:
					maxmissing_filtered = True
			if not maxmissing_filtered:
				diff['ref'] = ref
				diff['alts'] = alts
				fixed_diffs[(chrom, base)] = diff
				diffs_found +=1

# exit if input wasn't recognized
if not input_read:
	print("Could not recognize format of input file, must be a fasta (.fa,.fasta,.fa.gz,.fasta.gz) or a vcf (.vcf,.vcf.gz).")	

if diffs_found > 0:
	with open(args.output_prefix + ".tsv", 'w') as of:
		if not args.verbose_output:
			of.write('\t'.join(['pos'] + [pop for pop in pops.keys()]) + "\n")
			for pos in fixed_diffs.keys():
				of.write('\t'.join([str(pos[1])] + [fixed_diffs[pos][pop][0] for pop in pops.keys()]) + "\n")
			of.close()
		else:
			print("Writing verbose output file")
			of.write('\t'.join(['chrom','pos','ref','alt'] + [pop for pop in pops.keys()]) + "\n")
			for pos in fixed_diffs.keys():
				of.write(str(pos[0]) + "\t" + str(pos[1]) + "\t" + str(fixed_diffs[pos]['ref']) + "\t" + str(fixed_diffs[pos]['alts']) + "\t" + '\t'.join([fixed_diffs[pos][pop][0] for pop in pops.keys()]) + "\n")
			of.close()
	if args.output_vcf:
		print("Writing vcf output...")
		with open(args.output_prefix + ".vcf", 'w') as of:
			for pos in fixed_diffs.keys():
				of.write(str(pos[0]) + "\t" + str(pos[1]) + "\t" + "." + "\t" + fixed_diffs[pos]['ref'] + "\t" + fixed_diffs[pos]['alts'] + "\n")
	print("Wrote {count} fixed differences to {file}".format(count=diffs_found, file=args.output_prefix + ".tsv"))
else:
	print("No fixed differences were found between the specified populations - no output file will be written.")
