# Axel Jensen, 2022-10-21
import bgzip
import random
import pysam
# Some generally useful tools/classes for working with genomics data
def IUPAC(alleles):
	if alleles[0] == alleles[1] and alleles[0] in ['A','T','C','G','a','t','c','g']:
		return alleles[0]
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
		code = "N" * len(alleles[random.randint(0,1)])
	return code

# Class to store a sequence
class Sequence:
	def __init__(self,sample,sequence,metadata):
		self.sample = sample
		self.sequence = sequence
		self.length = len(sequence)
		self.metadata = metadata
	# function to convert a fasta string to a sequence object
	def fromFasta(self,fastastring, metadata=None):
		self.sample = fastastring.readlines()[0].rstrip().lstrip(">")
		self.sequence = ""
		for line in fastastring.readlines():
			if not line.startswith(">"):
				self.sequence += line.rstrip()
		self.length = len(self.sequence)
		self.metadata = metadata
	def addSequence(self,sequence):
		self.sequence = self.sequence + sequence
	def addIUPAC(self,hetcall):
		iupac = IUPAC(hetcall)
		self.addSequence(iupac)
	def addRandom(self,hetcall):
		self.addSequence(hetcall[random.randint(0,1)])


# class to store multiple sequences, either as an alignment or just a multi sequence object
class MultiSequence:
	def __init__(self, aligned=True):
		self.sequences = None
		self.alignment = aligned
		self.input_file = None
		self.output_file = None
		self.missing_intervals = None
		self.called_sites = None
	def vcfToMSA(self,vcf,chrom,start,end,samples,index=1,haploidize=True,heterozygotes="IUPAC",reference_genome=None,mask_missing_sites=True):
		if index == 1:
			start = start - 1
		records = vcf.fetch(chrom,start,end)
		# initiate a MultiSequence object including one sequence for each sample, or two if we're not haploidizing
		for sample in samples:
			if haploidize:
				self.sequences = [Sequence(sample,"",{'region':(chrom,start,end),'vcf_sample':sample}) for sample in samples]
				#msa = MultiSequence(sequences)
			else:
				haps1 = Sequence(sample + "__1","",{'region':(chrom,start,end),'haplotype':1,'vcf_sample':sample})
				haps2 = Sequence(sample + "__2","",{'region':(chrom,start,end),'haplotype':2,'vcf_sample':sample})
				self.sequences = haps1 + haps2
				#msa = MultiSequence(sequences)
		# keep track of the sites with called genotypes, for later projection onto the reference genome if wanted
		missing_intervals = []
		called_sites = []
		last_pos = start
		for rec in records:
			if rec.pos > last_pos + 1:
				missing_intervals.append((last_pos + 1,rec.pos - 1))
				called_sites.append(rec.pos)
				last_pos = rec.pos
			else:
				called_sites.append(rec.pos)
				last_pos = rec.pos
			alleles_dict = dict([(i,a) for i,a in enumerate(rec.alleles)]) # this gives a "translation" dict with integers as keys and alleles as values
			alleles_dict[None] = "N" # adds translation for missing genotypes too
			# if any alleles differ in length from each other (not recommended for now!), we make them the same length by adding dashes, and print a warning about this
			if not len(set([len(i) for i in alleles_dict.values()])) == 1:
				print("Warning: indel-type alleles of differing lengths found at position {}. This could lead to erroneous alignments.".format(rec.pos))
				longest_allele = 0
				for a in alleles_dict.values():
					if len(a) > longest_allele:
						longest_allele = len(a)
				for k in alleles_dict.keys():
					alleles_dict[k] = alleles_dict[k] + "-" * (longest_allele - len(alleles_dict[k]))
			# add the appropriate allele to each sequence object
			for seq in self.sequences:
				gt = rec.samples[seq.metadata['vcf_sample']]['GT']
				hetcall = (alleles_dict[gt[0]],alleles_dict[gt[1]])
				if haploidize:
					if None in gt:
						seq.addSequence("N")
						continue
					if heterozygotes == "IUPAC":
						seq.addIUPAC(hetcall)
						continue
					if heterozygotes == "random":
						seq.addRandom(hetcall)
				else:
					seq.addSequence(gt[seq.metadata['haplotype'] - 1])
		# if the last called position is smaller then end, it means we need to append another missing interval to the list
		if rec.pos < end:
			missing_intervals.append((rec.pos + 1,end))
		# if using a reference genome to extrapolate missing sites, do that now
		if reference_genome:
			# pull out the apropriate region from the refgenome
			refseq = reference_genome.fetch(chrom,start,end)
			# first make a little tuple list of missing vs called intervals
			intervals = []
			last_pos = start + 1
			for i in missing_intervals:
				if i[0] > last_pos:
					intervals.append((last_pos,i[0]-1,"called"))
					last_pos = i[1] +1
					intervals.append((i[0],i[1],"missing"))
				else:
					intervals.append((i[0],i[1],"missing"))
					last_pos = i[1] +1
			# if window ends with called bases, we need to append these to the list as well
			if intervals[-1][1] < end:
				intervals.append((intervals[-1][1] + 1, end, 'called'))
			# let's convert them to relative positions in this particular msa object, and 0-based
			intervals_con = []
			seqpos = 0
			for i in intervals:
				if i[2] == 'missing':
					relstart = i[0] - start - 1
					relend = i[1] - start
					intervals_con.append((relstart,relend,'missing'))
				else:
					relstart = seqpos
					relend = relstart
					for y in range(i[0] - 1, i[1]):
						relend += 1
						seqpos += 1
					intervals_con.append((relstart,relend,'called'))
		# now it should be straightforward to go through each sample, fetch the appropriate refsequence, and concatenate the called bases accordingly
		for seq in self.sequences:
			new_sequence = ""
			for i in intervals_con:
				if i[2] == 'missing':
					s = refseq[i[0]:i[1]]
					if mask_missing_sites:
						s = "N" * len(s)
				else:
					s = seq.sequence[i[0]:i[1]]
				new_sequence = new_sequence + s 
			seq.sequence = new_sequence
		#refseq = pysam.FastaFetch()
		self.missing_intervals = missing_intervals
		self.called_sites = called_sites
		#sprint(missing_intervals, called_sites)
	def length(self):
		# fetch all sequence lenghts
		lengths = [len(l.sequence) for l in self.sequences]
		if len(set(lengths)) == 1:
			return list(set(lengths))[0],"All sequences are of the same length."
		else:
			return lengths,"The sequences differ in length."
	def writePhylip(self,outfile,strict=False):
		with open(outfile, "w") as of:
			of.write('\t{nseq}\t{aln_length}\n'.format(nseq=len(self.sequences),aln_length=self.length()[0]))
			for seq in self.sequences:
				of.write(seq.sample + "   " + seq.sequence + "\n")
	def writeFasta(self,outfile,wrap=None):
		with open(outfile, "w") as of:
			for seq in self.sequences:
				of.write(">{}\n{}\n".format(seq.sample,seq.sequence))

# function that will initiate an msa object and extract add a vcf region to this
def vcfToMSA(vcf,chrom,start,end,samples=None,index=1,haploidize=True,heterozygotes="IUPAC",reference_genome=None,mask_missing_sites=True):
	if not samples:
		samples = list(vcf.header.samples)		
	# initiate an empty msa objec
	msa = MultiSequence()
	# pull out the alignment
	msa.vcfToMSA(vcf,chrom,start,end,samples,index=0,haploidize=haploidize,heterozygotes=heterozygotes,reference_genome=reference_genome,mask_missing_sites=mask_missing_sites)
	return msa