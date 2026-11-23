import gzip, sys
import array
vcf, chrom, parent, outfile = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
BIN = 2000000
with gzip.open(vcf, "rt") as fh:
    for line in fh:
        if line.startswith("#CHROM"):
            names = line.rstrip("\n").split("\t")[9:]
            iP = names.index(parent)
            n = len(names)
            match = {}
            total = {}
            break
    for line in fh:
        if not line.startswith(chrom + "\t"):
            if line[:6] > chrom:
                break
            continue
        cols = line.rstrip("\n").split("\t")
        pos = int(cols[1])
        if cols[4] not in ("A","C","G","T"):
            continue
        b = pos // BIN
        gts = [c.split(":")[0] for c in cols[9:]]
        pg = gts[iP].replace("|","/").split("/")
        if len(pg) != 2 or pg[0] == ".":
            continue
        if b not in match:
            match[b] = [0]*n; total[b] = [0]*n
        mb = match[b]; tb = total[b]
        ps = set(pg)
        for i, g in enumerate(gts):
            a = g.replace("|","/").split("/")
            if len(a) != 2 or a[0] == ".":
                continue
            tb[i] += 1
            if set(a) == ps:
                mb[i] += 1
with open(outfile, "w") as out:
    out.write("bin\tmedian_ibs\tmean_ibs\tn_high\tparent_ibs\n")
    for b in sorted(match):
        vals = []
        for i in range(n):
            if total[b][i] >= 50:
                vals.append(match[b][i] / total[b][i])
        if not vals:
            continue
        vals_sorted = sorted(vals)
        med = vals_sorted[len(vals)//2]
        mean = sum(vals)/len(vals)
        nhi = sum(1 for v in vals if v >= 0.9)
        pv = match[b][iP]/total[b][iP] if total[b][iP] else float("nan")
        out.write("%d\t%.4f\t%.4f\t%d\t%.4f\n" % (b, med, mean, nhi, pv))
print("bins:", len(match))
