# Axel Jensen 2022-10-25
# This script will resolve a gatk bug, where some ref alleles will span multiple nucleotide/positions erroneaously
from pysam import VariantFile
import argparse

def shorten_ref(rec):
	new_allele = rec.ref[0]
	new_rec = rec.copy()
	new_rec.alleles = (new_allele,".")
	return new_rec,rec

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",type=str,help="Input vcf to filter.")
parser.add_argument("-o","--output",type=str,help="Output vcf.")
parser.add_argument("--write-filtered-records",type=str,help="Optionally, give a path where to write a vcf with the records that's been filtered. Making it possible to restore the unfiltered vcf.", default=None)

args = parser.parse_args()

#input_vcf = "/crex/proj/sllstore2017021/nobackup/DENTI_LESULA_PROJECT/VARIANTS/MMUL_10/VCFs/FILTERED/ALLSITES/AUTOSOMES/allsites.filtered.NC_041772.1.mmul.vcf.gz"
#output_vcf = "./test_fix.vcf.gz"
#failed_vcf = "./test_failed.vcf.gz"
#write_filtered_records=True
input_vcf = args.input
output_vcf = args.output
failed_vcf = args.write_filtered_records

vcf = VariantFile(input_vcf)
outvcf = VariantFile(output_vcf,'w',header=vcf.header)
failvcf = VariantFile(failed_vcf,'w',header=vcf.header)

for rec in vcf.fetch():
	if len(rec.ref) == 1:
		outvcf.write(rec)
		continue
	elif len(rec.ref) > 1 and not rec.alts:
		new_rec,old_rec = shorten_ref(rec)
		outvcf.write(new_rec)
		if args.write_filtered_records:
			failvcf.write(old_rec)				
		else:
			outvcf.write(rec)
	elif len(rec.alts) > 0:
		outvcf.write(rec)

outvcf.close()
failvcf.close()
