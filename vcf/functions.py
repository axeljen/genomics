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

def getChromLengths(vcf, only_contigs_with_records=True):
	all_contigs = {contig: 0 for contig in list(vcf.header.contigs)}
	# fetch all the chromosomes to begin with from the vcf header
	for contig in all_contigs.keys():
		length = vcf.header.contigs[contig].length
		all_contigs[contig] = length
	contigs = {}
	if only_contigs_with_records:
		# usually we'll only want the chromosome lengths for the chromosome(s) that has records in this vcf file
		for contig,length in all_contigs.items():
			try: 
				# check this by trying to fetch the first record on each chromosome, if successfull, add this contig to the returned dict
				next(vcf.fetch(contig))
				contigs[contig] = length
			except:
				pass
	else:
		contigs = all_contigs 
	return contigs

#split a dictionary with genomic intervals into window of "window" size, separated by step_size
def generateWindows(intervals, window, step_size):
	# first, if the intervals are given as a simle chrom: length dictionary, add 0 as start position to each such interval
	if not type(intervals[list(intervals.keys())[0]]) is tuple:
		for c in intervals.keys():
			intervals[c] = (0, intervals[c])
    #create a dictionary with windows for each chromosome represented in intervals
	windows = []
	window_number = 1
	end = window
	for c in intervals.keys():
		start = intervals[c][0]
		length = intervals[c][1]
		chrom = c
		while start < length:
			windows.append({'window_number': window_number, 'chrom': chrom, 'start':start, 'end':end})
			start = start + step_size
			end = start + window if start + window < length else length
			window_number = window_number + 1
		else:
			break
	return windows

# function to parse a vcf file/region to a sequence dictionary
def vcfToSeqDict(vcf, chrom, haploidize = True, handle_heterozygotes = "IUPAC", samples=None, start=None, end=None):
	if not samples:
		samples = list(vcf.header.samples)
	if haploidize:
		seqDict = {s: [] for s in samples}
	else:
		seqDict = {}
		for s in samples:
			seqDict[s + "__1"] = []
			seqDict[s + "__2"] = []
	for rec in vcf.fetch(chrom,start,end):
		siteDict = {s: [] for s in seqDict.keys()}
		alleles = {0: rec.ref}
		c = 1
		for a in range(0,len(rec.alts)):
			alleles[c] = rec.alts[a]
			c +=1
		alleles['missing'] = "N"
		# add dashes to the alleles to match up the longest allele, in case it's an indel or mixed
		longest_allele = 0
		for a in alleles.values():
			if len(a) > longest_allele:
				longest_allele = len(a)
		if longest_allele > 1:
			print("Warning: found alleles longer then 1 bp (i.e. not only biallelic SNPs). This could lead to problems in the output alignment, as alleles of differing lengths may come out as non-aligned.")
		for k,a in alleles.items():
			allele = a + "-" * (longest_allele - 1 - len(a))
		for s in samples:
			gt_int = rec.samples[s]['GT']
			if None in gt_int:
				gt_dip = (alleles['missing'], alleles['missing'])
			else:
				gt_dip = (alleles[gt_int[0]], alleles[gt_int[1]])
			if haploidize:
				if len(set(gt_dip)) > 1:
					if handle_heterozygotes == "IUPAC":
						gt_hap = IUPAC(gt_dip)
					elif handle_heterozygotes == "random":
						gt_hap = gt_dip[random.randint(0,1)]
				else:
					gt_hap = gt_dip[0]
				siteDict[s].append(gt_hap)
			else:
				siteDict[s + "__1"] = [gt_dip[0]]
				siteDict[s + "__2"] = [gt_dip[1]]
		seqDict = {k: v + siteDict[k] for k,v in seqDict.items()}
	# last we join the list together as strings, and do a sanity check confirming all sequences are of the same lengths
	lengths = []
	for k,v in seqDict.items():
		seqDict[k] = ''.join(v)
		lengths.append(len(v))
	if len(set(lengths)) > 1:
		print("The sequences came out with differing lengths for some reason, which is not ok. Exiting.")
		sys.exit()
	return seqDict

def filterSeqDict(seqDict, maxmissing):
	lengths = [len(seqDict[k]) for k in seqDict.keys()]
	if len(set(lengths)) > 1:
		print("Differing lengths found in this alignment, exiting..")
		sys.exit()
	length = lengths[0]
	missing = {i: 0 for i in range(0, length)}
	for i in range(0,length):
		for s in seqDict.keys():
			if seqDict[s][i] in ['N','n','-']:
				missing[i] += 1
		missing[i] = missing[i] / len(list(seqDict.keys()))
	filtered_seqDict = {s: "" for s in seqDict.keys()}
	for pos,missing in missing.items():
		if missing <= maxmissing:
			for s in seqDict.keys():
				filtered_seqDict[s] = filtered_seqDict[s] + seqDict[s][pos]
	return filtered_seqDict

