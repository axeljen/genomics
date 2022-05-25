
def parsePopfile(popfile):
	clades = {}
	with open(popfile) as f:
		for line in f.readlines():
			sample = line.split("\t")[0]
			try:
				clade = line.split("\t")[1].rstrip()
				if not clade in list(clades.keys()):
					clades[clade] = [sample]
				else:
					clades[clade].append(sample)
			except:
				continue
	return clades

def DifferentiallyFixed(pops, seqDict):
	diff = False
	popalleles = {}
	for pop,samples in pops.items():
		popalleles[pop] = []
		for sample in samples:
			if not seqDict[sample][0] in popalleles[pop] and seqDict[sample][0] in ['A','T','C','G','a','t','c','g']:
				popalleles[pop].append(seqDict[sample][0])
	#check if the site is differentially fixed between at least two pops
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
		return popalleles
	else:
		return False

	
import seqtools.functions as seqtools
import argparse

parser = argparse.ArgumentParser()

#arguments for in and output:
parser.add_argument("-i", "--input-fasta", help="Input alignment in fasta format.", required=True)
parser.add_argument("-o", "--output-table", help="Output table in tsv format.", required=True)

parser.add_argument("-p", "--popfile", help="File with two tab separated columns, sample in the first and clade/population in second. It's ok to have samples without assignment.", required=True)

args = parser.parse_args()


#fasta = "/Users/axeljensen/Dropbox/FRANKIES_PROJECT/TFB2M.fa"

# read in fasta
seqDict = seqtools.FastaFetch(args.input_fasta)

# read popfile and parse it
popfile=args.popfile
pops = parsePopfile(popfile)

# check length of alignment
seqlen = len(seqDict[list(seqDict.keys())[0]])

# loop through each position to look for fixed differences
fixed_diffs = {}
diffs_found = 0
for base in range(0,seqlen):
	posDict = {sample: seqDict[sample][base] for sample in seqDict.keys()}
	diff = DifferentiallyFixed(pops, posDict)
	if diff:
		fixed_diffs[base + 1] = diff
		diffs_found +=1
with open(args.output_table, 'w') as of:
	of.write('\t'.join(['pos'] + [pop for pop in pops.keys()]) + "\n")
	for pos in fixed_diffs.keys():
		of.write('\t'.join([str(pos)] + [fixed_diffs[pos][pop][0] for pop in pops.keys()]) + "\n")
	of.close()

print("Wrote {count} fixed differences to {file}".format(count=diffs_found, file=args.output_table))