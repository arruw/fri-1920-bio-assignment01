import helpers as h

print("Fetching sequence...")
seq = h.getSequence()
print("Done.")

print("Fetching annotated genes...")
genes = h.getGenes()
print("Done.")

print("Finding possible proteins...")
prots = h.findProteins(1, seq) | h.findProteins(-1, h.reverseComplement(seq))
print(f"Done. Found {len(prots)} possible genes.")

Ps = list()
Rs = list()
Ls = range(50, 500, 5)

print("L\t| TP\tFP\t| P\tR\t| F1")
print("------------------------------------------------------")
for L in Ls:
    detected = set(filter(lambda x: (x[2] - x[1] + 1)/3 > L, prots))

    N = len(detected)
    TP = len(genes & detected)
    FP = N - TP
    P = TP / N
    R = TP / len(genes)
    F1 = 2*P*R/(P+R)

    Ps.append(P)
    Rs.append(R)

    print(f"{L}\t| {TP}\t{FP}\t| {P:.2f}\t{R:.2f}\t| {F1:.2f}")

import matplotlib.pyplot as plt

plt.plot(Ls, Ps)
plt.plot(Ls, Rs)
plt.xlabel("Length treshold")
plt.title("Precision & Recall compared to the length treshold")
plt.legend(['Precision', 'Recall'], loc='upper right')
plt.savefig("L-PR-plot.png")
plt.show()