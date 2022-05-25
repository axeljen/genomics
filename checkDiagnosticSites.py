

def parseSites(diagnostic_sites):
	diagnostic_dict = {}
	with open(diagnostic_sites) as f:
		header = [pop.rstrip() for pop in f.readline().split("\t")[1:]]
		for line in f.readlines():
			pos = line.split("\t")[0].rstrip()
			diagnostic_dict[pos] = {}
			for i in range(0, len(header)):
				diagnostic_dict[pos][header[i]] = line.split("\t")[i +1].rstrip()
	return diagnostic_dict

def attachSampleAlleles(diagnostic_dict, seqDict):
	pops = list(diagnostic_dict[list(diagnostic_dict.keys())[0]].keys())
	counts = {}
	for sample in seqDict.keys():
		counts[sample] = {pop: 0 for pop in pops}
		for pos in diagnostic_dict.keys():
			base = seqDict[sample][int(pos) - 1]
			for pop in counts[sample].keys():
				#print(pop)
				if base == diagnostic_dict[pos][pop]:
					counts[sample][pop] += 1
			diagnostic_dict[pos][sample] = base
	return diagnostic_dict,counts


import seqtools.functions as seqtools
import argparse

parser = argparse.ArgumentParser(description="Uses the output from findFixedDifferences.py to count the site patterns for individual samples in an alignment.")

# arguments for in and output:
parser.add_argument("-i", "--input-fasta", help="Input alignment in fasta format.", required=True)
parser.add_argument("-d", "--diagnostic-sites", help="Table with diagnostic sites, as output from findFixedDifferences.py.")

parser.add_argument("-o", "--output-prefix", help="Prefix for the two output files: prefix_nucleotides.tsv and prefix_counts.tsv", required=True)

args = parser.parse_args()

# read and parse input fasta
input_fasta=args.input_fasta
seqDict = seqtools.FastaFetch(input_fasta)

# read and parse diagnostic sites
diagnostic_sites=parseSites(args.diagnostic_sites)
print(diagnostic_sites)
# copy this to attach all the samples
diagnostic_copy = {k: v for k,v in diagnostic_sites.items()}
print(diagnostic_copy)
#print(diagnostic_sites)
# count sites for all sasmples
result_dict,counts = attachSampleAlleles(diagnostic_copy, seqDict)
print(diagnostic_sites)
#print(diagnostic_sites)
# write the results dict
with open(args.output_prefix + "_nucleotides.tsv", "w") as of:
	samples = [s for s in result_dict[list(result_dict.keys())[0]]]
	of.write('\t'.join(['pos'] + samples) + "\n")
	for site in result_dict.keys():
		of.write('\t'.join([site] + [result_dict[site][s] for s in samples]) + "\n")
	of.close()

# and then the counts
with open(args.output_prefix + "_counts.tsv", "w") as of:
	pops = [pop for pop in counts[list(counts.keys())[0]]]
	of.write('\t'.join(["sample"] + pops) + "\n")
	for sample in counts.keys():
		of.write('\t'.join([sample] + [str(counts[sample][pop]) for pop in pops]) + "\n")
	of.close()





#outprefix = "/Users/axeljensen/Dropbox/FRANKIES_PROJECT/diagnostic_test"


#diagnostic_dict = parseSites(diagnostic_sites)

#diagnostic_dict_e, counts = attachSampleAlleles(diagnostic_dict, seqDict)