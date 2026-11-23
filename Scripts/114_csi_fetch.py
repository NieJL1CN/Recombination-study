import gzip, struct, sys, zlib

def load_csi(path):
    with open(path, "rb") as f: head = f.read(2)
    data = gzip.open(path, "rb").read() if (head[0] == 31 and head[1] == 139) else open(path, "rb").read()
    assert data[:4] == b"CSI\x01"
    min_shift, depth, l_aux = struct.unpack_from("<iii", data, 4)
    off = 16 + l_aux
    n_ref = struct.unpack_from("<i", data, off)[0]; off += 4
    refs = []
    for _ in range(n_ref):
        n_bin = struct.unpack_from("<i", data, off)[0]; off += 4
        bins = {}
        for _ in range(n_bin):
            binid, loff, n_chunk = struct.unpack_from("<IQi", data, off); off += 16
            chunks = []
            for _ in range(n_chunk):
                u, v = struct.unpack_from("<QQ", data, off); off += 16
                chunks.append((u, v))
            bins[binid] = chunks
        refs.append(bins)
    aux = data[16:16 + l_aux]
    l_nm = struct.unpack_from("<i", aux, 28)[0]
    names = aux[32:32 + l_nm].split(b"\x00")
    return min_shift, depth, refs, [n.decode() for n in names if n]

def reg2bins(beg, end, min_shift, depth):
    end -= 1; out = [0]; t = 0; s = min_shift
    for l in range(depth + 1):
        b = t + (beg >> s); e = t + (end >> s)
        out.extend(range(b, e + 1)); t += 1 << (l * 3); s += 3
    return out

def bgzf_block(f, coffset):
    f.seek(coffset); head = f.read(12)
    if len(head) < 12: return None, None
    xlen = struct.unpack_from("<H", head, 10)[0]; extra = f.read(xlen)
    bsize = None; i = 0
    while i + 4 <= xlen:
        si1, si2, slen = struct.unpack_from("<BBH", extra, i)
        if si1 == 66 and si2 == 67:
            bsize = struct.unpack_from("<H", extra, i + 4)[0]; break
        i += 4 + slen
    rest = f.read(bsize + 1 - 12 - xlen)
    data = zlib.decompress((head + extra + rest)[18:], -15)
    return data, coffset + bsize + 1

def fetch(vcf, csi, chrom, beg, end, max_lines=500000):
    min_shift, depth, refs, names = load_csi(csi)
    rid = names.index(chrom)
    chunks = []
    for b in reg2bins(beg, end, min_shift, depth):
        chunks.extend(refs[rid].get(b, []))
    chunks.sort()
    out = []
    with open(vcf, "rb") as f:
        for (vstart, vend) in chunks:
            coff = vstart >> 16; uoff = vstart & 0xFFFF; coff_end = vend >> 16
            while True:
                data, nxt = bgzf_block(f, coff)
                if data is None: break
                seg = data[uoff:]
                for line in seg.split(b"\n"):
                    if not line.startswith(b"chr"): continue
                    cols = line.split(b"\t")
                    if len(cols) < 9: continue
                    if cols[0].decode() != chrom: continue
                    pos = int(cols[1])
                    if beg <= pos <= end:
                        out.append((pos, cols))
                        if len(out) >= max_lines: return out
                if nxt > coff_end: break
                coff = nxt; uoff = 0
    return out

if __name__ == "__main__":
    r = fetch(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
    print("records:", len(r))
    for pos, cols in r[:3]:
        print(pos, cols[3].decode(), cols[4].decode(), len(cols))