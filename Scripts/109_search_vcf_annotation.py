import gzip, sys
vcf = sys.argv[1]
targets = [t.encode() for t in sys.argv[2:]]
buf = b""
count = 0
hits = []
with gzip.open(vcf, "rb") as fh:
    while True:
        chunk = fh.read(1 << 23)
        if not chunk:
            break
        data = buf + chunk
        lines = data.split(b"\n")
        buf = lines[-1]
        for ln in lines[:-1]:
            for t in targets:
                if t in ln:
                    hits.append(ln[:400])
                    break
        if len(hits) > 40:
            break
print("scanned bytes ok; hits:", len(hits))
for h in hits[:40]:
    print(h.decode("utf-8", "ignore")[:300])
