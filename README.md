# genomics
Collection of python scripts used to process sequence data in various formats, and some basic statistics etc.

## Installation and dependencies

All the scripts run in python3, and require the dependencies listed in requirements.txt. Once you have python3 installed, Clone the repository and install the dependencies:

	git clone https://github.com/axeljen/genomics.git
	
	cd genomics

	python3 -m pip install -r requirements.txt


Then the scripts can be run directly out of this directory.

## VcfToMSA.py

	usage: VcfToMSA.py [-h] -i INPUT -o OUTPUT [-r REGION] [--samples SAMPLES] [--ignore-filters] [--set-filtered-to SET_FILTERED_TO] [--haploidize]
                   [--diploidize] [--handle-heterozygotes {iupac,random}] [-R REFERENCE_SEQUENCE] [--fill {reference,N}] [--verbose]

	This tool will read a VCF file and convert it to a multi sequence alignment.

	optional arguments:
	-h, --help            show this help message and exit
	-i INPUT, --input INPUT
							Input vcf
	-o OUTPUT, --output OUTPUT
							Output file in fasta (ending with .fa or .fa.gz) or in phylip (ending with .phy or phy.gz).
	-r REGION, --region REGION
							Region to extract as chrom:start-end, format is 1-based inclusive.
	--samples SAMPLES     File with samples to include, one per line. Defaults to all samples in vcf.
	--ignore-filters      If given, sites will be output regardless of filters in the INFO-column. Default is changing filtered sites to 'N'
	--set-filtered-to SET_FILTERED_TO
							Optionally give a custom letter for changing filtered sites to.
	--haploidize          Default behaviour: output a haploidized sequence for each sample.
	--diploidize          Will output two pseudo-haplotypes per sample, not taking phasing into account. Overrides '--haploidize'
	--handle-heterozygotes {iupac,random}
							Can be used in conjunction with '--haploidize': choices are 'iupac' (default) or 'random': random will sample a random allele in
							case of heterozygosity, iupac uses iupac ambiguity codes.
	-R REFERENCE_SEQUENCE, --reference-sequence REFERENCE_SEQUENCE
							If a reference genome is provided, the variants in the vcf will be output as an 'alignment' to this, with N's (default) or the
							reference bases filling out any gaps in the vcf.
	--fill {reference,N}  'reference' will use the reference bases to fill gaps not covered in the vcf, 'N' (default) will use Ns.
	--verbose             Will output some extra info like warnings and roadmarks for debugging.

Convert variants in a vcf file to a multi sequence alignment file in fasta or phylip format. You can choose between a few different behaviours:

* If a reference sequence in fasta format is supplied, then this can be used to fill gaps in between variants. This can be handy if working with bi-allelic SNPs for example, as it will not just concatenate the SNPs but also fill the gaps in between. The gaps can be filled either with the reference base or with N's, as chosen by the `--fill` flag.

* Default behavior is to merge diploid sequences into a single pseudo-haplotype (`--haploidized`). In this mode, heterozygote sites can be handled either by randomly sampling one allele or using IUPAC ambiguity codes, as specified with `--handle-heterozygotes`.

* Alternatively, you can choose to output two "haplotypes" per sample by specifying `--diploidize`. This option does not take phasing info into account at all but currently just puts the reference allele in one haplotype and the alternate in another (in the case of bi-allelic snps, otherwise it'll just take the first and second allele). So it's not very useful right now.

### Note: Use the fasta output format for now as the phylip is buggy. Will correct this soon..


## findFixedDifference.py

This script will search through an alignment and identify sites that are fixed for different alleles between specified groups.

	usage: findFixedDifferences.py [-h] -i INPUT [-r REGION] -o OUTPUT_TABLE -p POPFILE

	optional arguments:
	-h, --help            show this help message and exit
	-i INPUT, --input INPUT
							Input file, either as alignment in fasta format or as a vcf file.
	-r REGION, --region REGION
							Region to analyze as chrom:start-end, compatible only with vcf input.
	-o OUTPUT_TABLE, --output-table OUTPUT_TABLE
							Output table in tsv format.
	-p POPFILE, --popfile POPFILE
							File with two tab separated columns, sample in the first and clade/population in second. It's ok to have samples without assignment.

Takes a fasta sequence, or a vcf-file (or region) and a table specifying the populations/groups as input. The popfile should have one row per sample with two tab-separated columns: column 1 is sample name, and column 2 is population/group. It's ok to have samples unassigned. Output is a tab-separated table with one row for each position with a fixed difference, giving the alleles that each of the populations carry.

## checkDiagnosticSites.py

This tool does a follow up analysis on the findFixedDifferences.py run on a fasta alignment script by checking the site patterns in all samples in an alignment.

	usage: checkDiagnosticSites.py [-h] -i INPUT_FASTA [-d DIAGNOSTIC_SITES] -o OUTPUT_PREFIX

	Uses the output from findFixedDifferences.py to count the site patterns for individual samples in an alignment.

	optional arguments:
	-h, --help            show this help message and exit
	-i INPUT_FASTA, --input-fasta INPUT_FASTA
						Input alignment in fasta format.
	-d DIAGNOSTIC_SITES, --diagnostic-sites DIAGNOSTIC_SITES
						Table with diagnostic sites, as output from findFixedDifferences.py.
	-o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
						Prefix for the two output files: prefix_nucleotides.tsv and prefix_counts.tsv

