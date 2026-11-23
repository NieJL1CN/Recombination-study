import struct, zlib, sys, gzip

def read_csi(path):
    with open(path, 'rb') as f:
        head = f.read(2)
    if head[0] == 31 and head[1] == 139:
        with gzip.open(path, 'rb') as f:
            data = f.read()
    else:
        with open(path, 'rb') as f:
            data = f.read()
            data = f.read()
    assert data[:4] == b"CSI\x01", "not CSI"
    min_shift, depth, l_aux = struct.unpack_from("<iii", data, 4)
    off = 16
    aux = data[off:off+l_aux]; off += l_aux
    n_ref = struct.unpack_from("<i", data, off)[0]; off += 4
    refs = []
    for _ in range(n_ref):
        n_bin = struct.unpack_from("<i", data, off)[0]; off += 4
        bins = {}
        for _ in range(n_bin):
            b, loff, n_chunk = struct.unpack_from("<IQi", data, off); off += 16
            chunks = []
            for _ in range(n_chunk):
                u, v = struct.unpack_from("<QQ", data, off); off += 16
                chunks.append((u, v))
            bins[b] = chunks
        refs.append(bins)
    return min_shift, depth, aux, refs

def reg2bin(beg, end, min_shift, depth):
    end -= 1
    l = depth
    s = min_shift
    t = ((1 << (depth * 3)) - 1) // 7
    while l > 0:
        if beg >> s == end >> s:
            return t + (beg >> s)
        l -= 1; s += 3; t -= 1 << (l * 3)
    return 0

def reg2bins(beg, end, min_shift, depth):
    end -= 1
    out = [0]
    t = 0
    s = min_shift
    for l in range(depth + 1):
        b = t + (beg >> s); e = t + (end >> s)
        out.extend(range(b, e + 1))
        t += 1 << (l * 3)
        s += 3
    return out

def bgzf_read_block(f, coffset):
    f.seek(coffset)
    head = f.read(12)
    if len(head) < 12: return None, None
    xlen = struct.unpack_from("<H", head, 10)[0]
    extra = f.read(xlen)
    bsize = None
    i = 0
    while i + 4 <= xlen:
        si1, si2, slen = struct.unpack_from("<BBH", extra, i)
        if si1 == 66 and si2 == 67:
            bsize = struct.unpack_from("<H", extra, i+4)[0]
            break
        i += 4 + slen
    if bsize is None: raise RuntimeError("no BSIZE")
    rest = f.read(bsize + 1 - 12 - xlen)
    block = head + extra + rest
    data = zlib.decompress(block[18:], -15)
    return data, coffset + bsize + 1

def fetch(vcf, csi, chrom, beg, end, max_lines=200000):
    min_shift, depth, aux, refs = read_csi(csi)
    # tabix aux: format, col_seq, col_beg, col_end, meta, skip, l_nm, names
    fmt, col_seq, col_beg, col_end = struct.unpack_from("<iiii", aux, 0)
    l_nm = struct.unpack_from("<i", aux, 16)[0]
    names = aux[20:20+l_nm].split(b"\x00")
    names = [n.decode() for n in names if n]
    try:
        rid = names.index(chrom)
    except ValueError:
        return []
    bins = refs[rid]
    chunks = []
    for b in reg2bins(beg, end, min_shift, depth):
        if b in bins:
            chunks.extend(bins[b])
    if not chunks: return []
    chunks.sort()
    out = []
    with open(vcf, "rb") as f:
        for (vstart, vend) in chunks:
            coff = vstart >> 16; uoff = vstart & 0xFFFF
            coff_end = vend >> 16
            while coff <= coff_end:
                data, nxt = bgzf_read_block(f, coff)
                if data is None: break
                text = data[uoff:] if coff == (vstart >> 16) else data
                for line in text.split(b"\n"):
                    if not line: continue
                    try:
                        cols = line.split(b"\t")
                        if cols[0].decode() != chrom: continue
                        pos = int(cols[1])
                    except Exception:
                        continue
                    if beg <= pos <= end:
                        out.append((pos, cols))
                        if len(out) >= max_lines: return out
                coff = nxt
                uoff = 0
                if coff > coff_end: break
    return out

if __name__ == "__main__":
    vcf, csi, chrom, beg, end = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
    r = fetch(vcf, csi, chrom, beg, end)
    print("records:", len(r))
    for pos, cols in r[:3]:
        print(pos, cols[3].decode(), cols[4].decode())