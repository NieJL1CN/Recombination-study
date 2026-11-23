import gzip, sys
vcf = sys.argv[1]
with gzip.open(vcf, "rt") as fh:
    for line in fh:
        if line.startswith("#CHROM"):
            names = line.rstrip("\n").split("\t")[9:]
            print("n samples:", len(names))
            for t in ["paragon","WATDE0019","WATDE0075","Chinese_Spring"]:
                print(t, "index:", names.index(t) if t in names else "MISSING")
            break
