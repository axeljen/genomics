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
