import gzip, sys
vcf, chrom, beg, end, parent, outfile = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6]
match = None; total = None; names = None
with gzip.open(vcf, "rt") as fh:
    for line in fh:
        if line.startswith("#CHROM"):
            names = line.rstrip("\n").split("\t")[9:]
            iP = names.index(parent)
            match = [0] * len(names); total = [0] * len(names)
            break
    for line in fh:
        if not line.startswith(chrom + "\t"):
            if line[:5] > chrom:
                break
            continue
        cols = line.rstrip("\n").split("\t")
        pos = int(cols[1])
        if pos < beg: continue
        if pos > end: break
        if cols[4] not in ("A","C","G","T"): continue
        gts = [c.split(":")[0] for c in cols[9:]]
        pg = gts[iP]
        for i, g in enumerate(gts):
            if g == "./." or "|" in g and False: 
                continue
            a = g.replace("|", "/").split("/")
            if len(a) != 2 or a[0] == ".": continue
            pa = pg.replace("|", "/").split("/")
            if len(pa) != 2 or pa[0] == ".": continue
            total[i] += 1
            if set(a) == set(pa): match[i] += 1
with open(outfile, "w") as out:
    out.write("sample\tibs\tn_sites\n")
    for i, n in enumerate(names):
        if total[i] > 50:
            out.write("%s\t%.4f\t%d\n" % (n, match[i] / total[i], total[i]))
print("done", len(names))
