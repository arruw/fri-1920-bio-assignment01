import helpers as h

seq = h.getSequence()
genes = h.getGenes()
# seq = 'ATG---TAG---ATG---ATG---+++TAA=ATG---TAG---ATG---ATG---+++TAA==ATG---TAG---ATG---ATG---+++TAA'


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
                proteins.add((strand, offset+start+1, offset+i+3))
                break

        # remove processed chunks
        rna = rna[start+3:]
        offset += start+3

    return proteins

prots = findProteins(1, seq) | findProteins(-1, reverseComplement(seq))

# print(prots)

print("L\t| P\tR\t| F1")
print("--------------------------------")
for L in range(50, 500, 5):
    detected = set(filter(lambda x: (x[2] - x[1] + 1)/3 > L, prots))

    nTP = len(genes & detected)
    P = nTP / len(detected)
    R = nTP / len(genes)
    F1 = 2*P*R/(P+R)

    print(f"{L}\t| {P:.2f}\t{R:.2f}\t| {F1:.2f}")