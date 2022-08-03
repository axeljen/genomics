# windowPhylogenies.py

This is a python script that can be used to construct phylogenies in sliding windows along the genome, taking a vcf file as input. It is also possible to provide a file with custom regions. For now, there's code to run iqtree and phyml, and these (or at least the one you want to use) need to installed. If it is not in your path, you can give the full path to the executable as an argument. It's using functions from sibling directories in this repo, so you'll need to clone the whole repo to be able to run this.

	usage: windowPhylogenies.py [-h] -i snps.vcf | spns.vcf.gz -o output.tsv [--sequences <seq1,seq2,seq3>|seqs.file]
                            [--max-missing <int:nBases>] [-w bases] [--step-size bases] [--haploidize]
                            [--handle-heterozygotes HANDLE_HETEROZYGOTES] [--iqtree-path path/to/iqtree] [--model MODEL]
                            [--outgroup outgroup_seq1,outgroup_seq2] [--bootstraps number_of_bootstraps] [-T threads]
                            [-R REGIONS] [--phyml] [--phyml-path PHYML_PATH] [--phyml-model PHYML_MODEL]
                            [--optimize OPTIMIZE] [--analysis-directory ANALYSIS_DIRECTORY] [--test]