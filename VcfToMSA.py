import argparse
import sys
import seqtools.functions as seqtools
from seqtools.functions import FastaFetch
import pysam
import vcf.functions as v
import os


parser = argparse.ArgumentParser(description="This tool will read a VCF file and convert it to a multi sequence alignment.")

# input and output
parser.add_argument('-i', '--input', type=str, help='Input vcf', required=True)
parser.add_argument('-o', '--output', type=str, help='Output file in fasta (ending with .fa or .fa.gz) or in phylip (ending with .phy or phy.gz).', required=True)

# regions
parser.add_argument('-r', '--region', type=str, help='Region to extract as chrom:start-end, format is 1-based inclusive.')

# samples
parser.add_argument('--samples', type=str, help="File with samples to include, one per line. Defaults to all samples in vcf.")

# handle filters
parser.add_argument('--ignore-filters', action="store_true",help="If given, sites will be output regardless of filters in the INFO-column. Default is changing filtered sites to 'N'")
parser.add_argument('--set-filtered-to', type=str, default="N", help="Optionally give a custom letter for changing filtered sites to.")

# diploidize or haploidize?
parser.add_argument('--haploidize', action="store_true", default=True, help="Default behaviour: output a haploidized sequence for each sample.")
parser.add_argument('--diploidize', action="store_true", help="Will output two pseudo-haplotypes per sample, not taking phasing into account. Overrides '--haploidize'")
parser.add_argument('--handle-heterozygotes', default="iupac", choices=("iupac","random"), help="Can be used in conjunction with '--haploidize': choices are 'iupac' (default) or 'random': random will sample a random allele in case of heterozygosity, iupac uses iupac ambiguity codes.")

# use reference genome?
parser.add_argument('-R', '--reference-sequence', type=str, help="If a reference genome is provided, the variants in the vcf will be output as an 'alignment' to this, with N's (default) or the reference bases filling out any gaps in the vcf.")
parser.add_argument('--fill', choices=("reference","N"), help="'reference' will use the reference bases to fill gaps not covered in the vcf, 'N' (default) will use Ns.")


parser.add_argument('--verbose', action="store_true", help="Will output some extra info like warnings and roadmarks for debugging.")
args = parser.parse_args()

# parse command inputs
input_vcf = args.input
outfile = args.output
if args.region:
	chrom,interval=args.region.split(":")[0],args.region.split(":")[1]
	start=int(interval.split("-")[0]) - 1
	end=int(interval.split("-")[1])
	region = [chrom,start,end]
else:
	region = [None,None,None]
if args.samples:
	samples = []
	with open(args.samples) as f:
		for line in f.readlines():
			if not line == "":
				samples.append(line.rstrip())
else:
	vcf = pysam.VariantFile(input_vcf)
	samples = list(vcf.header.samples)
if args.ignore_filters:
	ignore_filters = True
else:
	ignore_filters = False
set_filtered_to = args.set_filtered_to
if not args.diploidize:
	diploidize = False
	if args.handle_heterozygotes == "iupac":
		use_ambiguities = True
	elif args.handle_heterozygotes == "random":
		use_ambiguities = False
	else:
		sys.stderr.write("Error: invalid value given to --handle-heterozygotes, exiting...")
		sys.exit(1)

else:
	diploidize = True

# assume a reference is not given
reference = False

# if it is, read this into a seqdict
if args.reference_sequence:
	reference = True
	fill = args.fill
	refseq = FastaFetch(args.reference_sequence, sequences=chrom, start=start, end=end)

# Now we can start building the alignment
sys.stderr.write("Starting to process variants.\n")


# If we're using reference bases to fill the gaps, we start with a seqdict of the reference sequence for all samples, and replace it with the variants as we go
if reference:
	seqDict = {}
	refname = list(refseq.keys())[0]
	seqDict[refname] = refseq[refname]
	for sample in samples:
		if fill == "reference":
			if diploidize:
				seqDict[sample + "__1"] =  seqDict[refname]
				seqDict[sample + "__2"] = seqDict[refname]
			else:
				seqDict[sample] = seqDict[refname]
		else:
			if diploidize:
				seqDict[sample + "__1"] =  "N" * len(seqDict[refname])
				seqDict[sample + "__2"] = "N" * len(seqDict[refname])
			else:
				seqDict[sample] = "N" * len(seqDict[refname])
	# Keeping the sequences as a list at this point, for easier manipulation downstream
	for key,value in seqDict.items():
		seqDict[key] = [i for i in value]
else:
	seqDict = {sample: [] for sample in samples}
	seqDict['reference'] = []

# If we're filling out gaps with reference bases, we need to keep track of indels to get the offset position
offset = 0

# Keep track of how many variants we've processed to give some updates
processed_positions = 0
applied_variants = 0

with pysam.VariantFile(input_vcf) as vcf:
	for rec in vcf.fetch(region[0],region[1],region[2]):
		# first get all the alleles
		bases = v.get_bases(rec, samples, set_filtered_to=set_filtered_to, ignore_filters=ignore_filters)
		# and turn those in to a single-locus seqDict
		if diploidize:
			siteSeqDict = v.Diploidize(bases)
		else:
			siteSeqDict = v.Haploidize(bases)
		# if we're filling the gaps with reference or Ns, we need to keep track of the position to swap
		if reference:
			pos = rec.pos
			allele_length = max(len(a) for a in siteSeqDict.values())
			offset = offset + allele_length - 1
			seqDictPos = pos - start - 1 + offset
			if allele_length > 1: # if there's an insertion we need to adjust the reference sequence with gaps
				refallele = seqDict[refname][seqDictPos] + "-" * (allele_length - len(seqDict[refname][seqDictPos]))
				seqDict[refname][seqDictPos] = refallele
				for sample in siteSeqDict.keys():
					seqDict[sample][seqDictPos] = siteSeqDict[sample] + "-" * (allele_length - len(siteSeqDict[sample]))
			else: #if everything is of length 1 just replace the correct position with the new base
				for sample in siteSeqDict.keys():
					seqDict[sample][seqDictPos] = siteSeqDict[sample]
		# And if we don't care about filling out the gaps (i.e. only concatenating variants), it's a simpler procedure:
		else:
			seqDict['reference'].append(rec.alleles[0])
			for sample in samples:
				seqDict[sample].append(siteSeqDict[sample])
			seqDictPos = len(seqDict['reference'])
		applied_variants += 1
		processed_positions = seqDictPos

# convert the sequences to strings
for key,value in seqDict.items():
	seqDict[key] = ''.join(value)

# write seqdict to output file
sys.stderr.write("Writing output alignment to {}\n".format(outfile))

# determine the output format
phylip = outfile.endswith((".phy", ".phylip", ".phy.gz", ".phylip.gz"))
fasta = outfile.endswith((".fa",".fa.gz",".fasta",".fasta.gz"))

if phylip:
	seqtools.WriteSeqDictToPhylip(seqDict, outfile, strict=False, append=False)
elif fasta:
	seqtools.WriteSeqDictToFasta(seqDict, outfile)

# write summary and exit
sys.stderr.write("Finished. Applied {variants} over a total of {seqlen} bases.\n".format(variants = applied_variants, seqlen = processed_positions))
