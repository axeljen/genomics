import pysam, os
from Bio import bgzf

#Function to fo from a diploid genotype (two alleles) to a IUPAC code
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

#function to parse a fastafile to seqdict, requires pysam FastaFile to be loaded
def FastaFetch(fastafile, sequences=None, start=None, end=None):
    f = pysam.FastaFile(fastafile)
    seqDict = {}
    if sequences is not None:
        if os.path.isfile(sequences):
            seqList = []
            with open(sequences) as s:
                for line in s.readlines():
                    if not line == "":
                        seqList.append(line.rstrip().lstrip(">")) 
        else:
            seqList = []
            for s in sequences.split(","):
                if not s == "":
                    seqList.append(s.rstrip().lstrip(">"))
        if start and end:
            for seq in seqList:
                seqDict[seq] = f.fetch(seq, start, end)
        elif start:
            for seq in seqList:
                seqDict[seq] = f.fetch(seq, start)
        else:
            for seq in seqList:
                seqDict[seq] = f.fetch(seq)
    else:
        if start and end:
            seqDict = []
            for seq in f.references:
                seqDict[seq] = f.fetch(seq, start, end)
        elif start:
            for seq in f.references:
                seqDict[seq] = f.fetch(seq, start)
        else:
            for seq in f.references:
                seqDict[seq] = f.fetch(seq)
    return seqDict

def openFile(filepath, mode='r', append=False):
    if filepath.endswith(".gz"):
        if append:
            return bgzf.open(filepath, 'rb')
        else:
            return bgzf.open(filepath, 'wt')
    else:
        if append:
            return open(filepath, 'a')
        else:
            return open(filepath, 'w')

def closeFile(f, filepath):
    if filepath.endswith(".gz"):
        return bgzf.close(f)
    else:
        return f.close()

def WriteSeqDictToFasta(seqDict, outfile, append=False):
    with openFile(outfile) as f:
        for name,seq in seqDict.items():
            f.write('>' + name + "\n" + seq + "\n")
        closeFile(f, outfile)

def WriteSeqDictToPhylip(seqDict, outfile, strict=False, append=False):
    with openFile(outfile, append=append) as f:
        f.write('\t' + str(len(seqDict)) + "\t" + str(len(seqDict[list(seqDict.keys())[0]])) + "\n")
        for key,value in seqDict.items():
            if not strict:
                header = key
            else:
                if len(key) < 10:
                    header = key + str(" " * 10 - len(key))
                else:
                    header = key[0:10]
            f.write(header + "   " + value + "\n")
        f.close()