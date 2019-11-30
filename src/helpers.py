__EMAIL = "mm3058@student.uni-lj.si"
__ID = "NC_000908"

def getSequence():
    from Bio import Entrez, SeqIO
    Entrez.email = __EMAIL

    with Entrez.efetch(db="nucleotide", id=__ID,rettype="fasta") as h:
        seq = str(SeqIO.read(h, "fasta").seq)
        return seq

def getGenes():
    from Bio import Entrez, SeqIO
    from statistics import median
    Entrez.email = __EMAIL

    with Entrez.efetch(db="nucleotide", id=__ID, rettype="gb") as h:
        gb = SeqIO.read(h, "gb")
        genes = list(filter(lambda f: f.type == "CDS" and 'translation' in f.qualifiers, gb.features))

        lengths = list(map(lambda f: (f.location.nofuzzy_end - f.location.nofuzzy_start)/3, genes))

        print(f"Shortest length: {min(lengths)}")
        print(f"Longest length:  {max(lengths)}")
        print(f"Median length:   {median(lengths)}")

        sst = set(map(lambda f: (f.strand, f.location.nofuzzy_start+1, f.location.nofuzzy_end), genes))
        return sst

def reverseComplement(rna):
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
    }
    return ''.join([complement[c] for c in list(rna[::-1])])

def findProteins(strand, rna):
    stop = {
        'TAA': True,
        'TAG': True
    }

    offset = 0
    proteins = set()
    while(True):
        # find start codon
        start = rna.find('ATG')
        if start == -1: break

        for i in range(start, len(rna), 3):
            # extract codon
            codon = rna[i:i+3]

            # we are at the end of sequence
            if len(codon) != 3: break

            # we found the stop codon
            if stop.get(codon, False):
                start_i = offset+start+1
                end_i = offset+i+3
                if strand == 1:
                    proteins.add((strand, start_i, end_i))
                elif strand == -1:
                    proteins.add((strand, 580076-end_i+1, 580076-start_i+1))
                break

        # remove processed chunks
        rna = rna[start+3:]
        offset += start+3

    return proteins

