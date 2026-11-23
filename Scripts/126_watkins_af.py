import gzip, sys
vcf, chrom, beg, end, outfile = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
with gzip.open(vcf, "rt") as fh:
    for line in fh:
        if line.startswith("#CHROM"):
            names = line.rstrip("\n").split("\t")[9:]
            wat = [i for i, n in enumerate(names) if n.startswith("WATDE")]
            others = {n: i for i, n in enumerate(names) if n in ("paragon", "Chinese Spring", "Chinese_Spring")}
            ip = others.get("paragon", others.get("Chinese_Spring", -1))
            break
    out = open(outfile, "w", encoding="utf-8")
    out.write("chr\tpos\tref\talt\tn_wat\twat_af\tn_hom_ref\tn_het\tn_hom_alt\tn_miss\n")
    n = 0
    for line in fh:
        if not line.startswith(chrom + "\t"):
            continue
        cols = line.rstrip("\n").split("\t")
        pos = int(cols[1])
        if pos < beg or pos > end:
            continue
        gts = [c.split(":")[0] for c in cols[9:]]
        af = 0.0; nr = nh = na = nm = 0
        for i in wat:
            g = gts[i].replace("|", "/")
            if g in ("./.", ".", ""):
                nm += 1; continue
            a = g.split("/")
            if len(a) != 2:
                nm += 1; continue
            s = sum(1 for x in a if x == "1")
            af += s / 2.0
            if s == 0: nr += 1
            elif s == 1: nh += 1
            else: na += 1
        ncall = nr + nh + na
        if ncall == 0:
            continue
        out.write("%s\t%d\t%s\t%s\t%d\t%.4f\t%d\t%d\t%d\t%d\n" % (
            cols[0], pos, cols[3], cols[4], ncall, af / ncall, nr, nh, na, nm))
        n += 1
    out.close()
print("sites:", n, "watkins samples:", len(wat))
