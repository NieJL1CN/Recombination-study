import gzip, sys
vcf, chrom, beg, end, outfile = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
samples = ["paragon", "robigus", "WATDE0019", "Chinese_Spring"]
with gzip.open(vcf, "rt") as fh:
    for line in fh:
        if line.startswith("#CHROM"):
            hdr = line.rstrip("\n").split("\t")
            idx = [hdr.index(s) for s in samples]
            break
    out = open(outfile, "w", encoding="utf-8")
    out.write("chr\tpos\tref\talt\tgt_paragon\tgt_robigus\tgt_W0019\tgt_CS\n")
    n = 0
    for line in fh:
        if not line.startswith(chrom + "\t"):
            continue
        cols = line.rstrip("\n").split("\t")
        pos = int(cols[1])
        if pos < beg or pos > end:
            continue
        gts = [cols[i].split(":")[0] for i in idx]
        out.write("\t".join([cols[0], str(pos), cols[3], cols[4]] + gts) + "\n")
        n += 1
    out.close()
print("variants:", n)
