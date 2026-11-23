import gzip, sys
vcf = sys.argv[1]
targets = [163839715, 163840445, 163840600, 163844000]
with gzip.open(vcf, "rt") as fh:
    hdr = None
    for line in fh:
        if line.startswith("#CHROM"):
            hdr = line.rstrip("\n").split("\t"); break
    names = hdr[9:]
    ip = names.index("paragon"); i7 = names.index("WATDE0075"); ics = names.index("Chinese_Spring")
    print("indices:", ip, i7, ics, "| sample check:", names[ip], names[i7], names[ics])
    hits = 0
    for line in fh:
        if not line.startswith("chr4B\t"): continue
        cols = line.rstrip("\n").split("\t")
        pos = int(cols[1])
        if pos not in targets: continue
        gts = [c.split(":")[0] for c in cols[9:]]
        alt = sum(1 for g in gts if "1" in g)
        print(pos, cols[3], cols[4], "| paragon", gts[ip], "| W579", gts[i7], "| CS", gts[ics],
              "| ALT carriers:", alt, "/", len(gts))
        hits += 1
        if hits >= len(targets): break
