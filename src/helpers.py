__EMAIL = "mm3058@student.uni-lj.si"
__ID = "NC_000908"

def getSequence():
    try:
        with open("store/sequence.txt", "r") as f:
            return f.read()
    except:
        from Bio import Entrez, SeqIO
        Entrez.email = __EMAIL

        with Entrez.efetch(db="nucleotide", id=__ID,rettype="fasta") as h:
            seq = str(SeqIO.read(h, "fasta").seq)
            with open("sequence.txt", "w") as f:
                return f.write(seq)
            return seq

def getGenes():
    try:
        with open("store/genes.txt", "r") as f:
            return f.read()
    except:
        from Bio import Entrez, SeqIO
        from statistics import median
        Entrez.email = __EMAIL

        with Entrez.efetch(db="nucleotide", id=__ID, rettype="gb") as h:
            gb = SeqIO.read(h, "gb")
            
            # print("Size of genom: ", gb.features[0].location.nofuzzy_end/3)

            genes = list(filter(lambda f: f.type == "gene", gb.features))

            # print("Number of genes: ", len(list(genes)))

            # geneLen = list(map(lambda f: (f.location.nofuzzy_end-f.location.nofuzzy_start)/3, genes))

            # print("Shortest gene: ", min(geneLen))
            # print("Longest gene: ", max(geneLen))
            # print("Median length of gene: ", median(geneLen))


            sst = set(map(lambda f: (f.strand, f.location.nofuzzy_start+1, f.location.nofuzzy_end), genes))

            return sst