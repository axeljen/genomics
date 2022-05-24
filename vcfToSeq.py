
import argparse, sys

from pysam import VariantFile, FastaFile

from seqtools.functions import IUPAC, FastaFetch

parser = argparse.ArgumentParser()

#arguments for in and output:
parser.add_argument("-i", "--input-vcf", help="Input vcf.", required=True, metavar="variants.vcf | variants.vcf.gz")
parser.add_argument("-o", "--output", help="Output file.", required=True, metavar="output.fa|output.fa.gz")

parser.add_argument("--region", help="Genomic regions to include, given as chr:start-stop.", metavar="chrom:start-stop")
parser.add_argument("--samples", help="samples to include.", metavar="<seq1,seq2,seq3>", required=False)
parser.add_argument("--non-var", action="store_true", help="Include nonvarsites in alignment as 'N'")

args = parser.parse_args()


