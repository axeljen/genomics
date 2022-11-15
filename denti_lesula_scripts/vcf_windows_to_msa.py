import argparse
import sys
import os
dir_path = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, dir_path)
import genomics
import pysam
#set up and parse input arguments

parser = argparse.ArgumentParser(description="Extract and concatenate a list of regions from a vcf file into a multisequence alignment.")

parser.add_argument('--vcfs', type=str, help='Two col list of vcf files, first with path and second which chrom its connected to')

parser.add_argument('--out', type=str, help='Path to write concatenated alignment too.')

parser.add_argument('--windows', type=str, help="File listing the windows to concatenate in the correct order in three columns: chrom, start, end. Expecting 1-based inclusive indexed values.")

parser.add_argument('--samples', type=str, help="File with one sample per line, for all samples to include in the extraction.", default=None)

args = parser.parse_args()

vcf_paths = {line.split("\t")[1].rstrip(): line.split("\t")[0] for line in open(args.vcfs).readlines()}
windows = [{'chrom':line.split("\t")[0],'start':int(line.split("\t")[1]),'end':int(line.split("\t")[2].rstrip())} for line in open(args.windows).readlines()]
# add vcf path to each window as well
for w in windows:
	w['vcf'] = vcf_paths[w['chrom']]
# initiate an MSA object from the first region, and subsequently we'll add all windows to this object
sys.stderr.write("Concatenating {} windows...\n".format(len(windows)))
vcf = pysam.VariantFile(windows[0]['vcf'])
chrom,start,end = windows[0]['chrom'],windows[0]['start'],windows[0]['end']
# check if we're specifying samples:
samples = [s.rstrip() for s in open(args.samples).readlines() if not s == ""] if args.samples else None

msa = genomics.vcfToMSA(vcf,chrom,start,end,samples=samples,index=1,haploidize=True,heterozygotes="IUPAC",reference_genome=None,mask_missing_sites=True)

# then loop through the rest of them, and add to this object
counter = 0
for w in windows[1:]:
	vcf = pysam.VariantFile(w['vcf'])
	chrom,start,end = w['chrom'],w['start'],w['end']
	window_msa = genomics.vcfToMSA(vcf,chrom,start,end,samples=samples,index=1,haploidize=True,heterozygotes="IUPAC",reference_genome=None,mask_missing_sites=True)
	for seq in msa.sequences:
		wseq = window_msa.getSequenceByName(seq.sample)
		seq.addSequence(wseq.sequence)
	counter += 1
	if counter % 1 == 0:
		sys.stderr.write("{} windows concatenated.\n".format(counter))

# write the concatenated alignment
if args.out.endswith(("phy","phylip","phy.gz","phylip.gz")):
	msa.writePhylip(args.out)
elif args.out.endswith(("fa",".fasta",".fa.gz",".fasta.gz")):
	msa.writeFasta(args.out)
sys.stderr.write("Done. Concatenated alignment written to {}.\n".format(args.out))

	

