import gzip, sys
vcf = sys.argv[1]
with gzip.open(vcf, "rt") as fh:
    for line in fh:
        if line.startswith("#CHROM"):
            names = line.rstrip("\n").split("\t")[9:]
            hits = [n for n in names if "robig" in n.lower()]
            print("robigus-like:", hits)
            for t in ["paragon","WATDE0019","Chinese_Spring"]:
                print(t, names.index(t) if t in names else "MISSING")
            break
