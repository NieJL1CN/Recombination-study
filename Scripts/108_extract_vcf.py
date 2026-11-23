import gzip, sys, io
vcf, chrom, winfile, outfile = sys.argv[1:5]
wins = []
with open(winfile) as fh:
    next(fh)
    for line in fh:
        g, s, e = line.rstrip("\n").split("\t")
        wins.append((g, int(s), int(e)))
targets = ["paragon", "WATDE0075", "WATDE0100", "Chinese_Spring"]
idx = []
with gzip.open(vcf, "rt") as fh:
    for line in fh:
        if line.startswith("#CHROM"):
            hdr = line.rstrip("\n").split("\t")
            idx = [hdr.index(t) for t in targets]
            break
    out = io.open(outfile, "w", encoding="utf-8")
    out.write("gene\tchr\tpos\tref\talt\tgt_paragon\tgt_W579\tgt_W729\tgt_CS\n")
    n = 0
    for line in fh:
        if not line.startswith(chrom + "\t"):
            continue
        cols = line.rstrip("\n").split("\t")
        pos = int(cols[1])
        gene = None
        for g, s, e in wins:
            if s <= pos <= e:
                gene = g; break
        if gene is None:
            continue
        gts = [cols[i].split(":")[0] for i in idx]
        out.write("\t".join([gene, cols[0], str(pos), cols[3], cols[4]] + gts) + "\n")
        n += 1
    out.close()
print("variants:", n)
