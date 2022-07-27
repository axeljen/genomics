import random

# Function to go from a diploid genotype (two alleles) to a IUPAC code
def IUPAC(alleles):
	IUPAC = {
		'M': ['A','C'],
		'R': ['A','G'],
		'W': ['A','T'],
		'S': ['C','G'],
		'Y': ['C','T'],
		'K': ['G','T']
	}
	hit = False
	for c in IUPAC:
		if set(IUPAC[c]) == set(alleles):
			code = c
			hit = True
	if hit == False:
		code = "N"
	return code

# opposite direction, going from a iupac code to a tuple of bases
def IUPACReverse(allele):
	IUPAC = {
		'M': ['A','C'],
		'R': ['A','G'],
		'W': ['A','T'],
		'S': ['C','G'],
		'Y': ['C','T'],
		'K': ['G','T']
	}
	hit = False
	for code in IUPAC.keys():
		if code == allele:
			bases = tuple(IUPAC[code])
	return bases
	



# this function will read a vcf record and return a tuple of each samples alleles
def get_bases(rec, samples=None, set_filtered_to="N", ignore_filters=False):
	sample_alleles = {}
	if samples==None:
		samples = list(rec.samples)
	alleles = {i: allele for i,allele in enumerate(rec.alleles)}
	longest_allele = 1
	for a in alleles:
		#in a nonvariant record, a lingering <NON_REF> might turn up at few places. These are low confident non-ref and should be removed, will set this to "N"
		if alleles[a] == "<NON_REF>":
			alleles[a] = "N"
	for a in alleles:
		if len(alleles[a]) > longest_allele:
			longest_allele = len(alleles[a])
		if len(alleles[a]) < longest_allele:
			alleles[a] = alleles[a] + "-" * (longest_allele - len(alleles[a]))
		if alleles[a] == "*" or alleles[a] == None:
			alleles[a] = "-" * longest_allele
	if not ''.join(list(rec.filter)) == 'PASS':
		if not ignore_filters:
			for sample in samples:
				sample_alleles[sample] = (set_filtered_to,set_filtered_to)
			return sample_alleles
	for sample in samples:
		adip = rec.samples[sample]['GT']
		if None in adip:
			a = "-" * longest_allele
			sample_alleles[sample] = (a,a)
		else:
			sample_alleles[sample] = (alleles[adip[0]], alleles[adip[1]])
	return sample_alleles

# this function will go through the bases of each sample (format {'sample':(allele1,allele2)}) and make a haploidized sequence from it
def Haploidize(sample_alleles, use_ambiguities=True):
	for sample in sample_alleles.keys():
		if use_ambiguities:
			varying_length = False
			for i in sample_alleles[sample]:
				last_len = len(i)
				if not len(i) == last_len:
					varying_length = True
			if varying_length:
				if args.verbose:
					print("Warning: " + sample + " is heterozygote for variants of varying lengths at " + str(rec.pos) + ", will randomly choose one of the alleles to output. ")
				a = sample_alleles[sample][random.randint(0,1)]
			else:
				if sample_alleles[sample][0] != sample_alleles[sample][1]:
					a = IUPAC(sample_alleles[sample])
				else:
					a = sample_alleles[sample][0]
			sample_alleles[sample] = a
		else:
			a = sample_alleles[sample][random.randint(0,1)]
			sample_alleles[sample] = a
	return sample_alleles
	
# Make a diploid sequence dictionary out of the vcf record alleles
def Diploidize(sample_alleles):
	seqDict = {}
	for sample in sample_alleles.keys():
		seqDict[sample + "__1"] = sample_alleles[sample][0]
		seqDict[sample + "__2"] = sample_alleles[sample][1]
	return seqDict

