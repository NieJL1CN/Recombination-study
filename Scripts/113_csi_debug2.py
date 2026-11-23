import gzip, struct, sys
path = sys.argv[1]
with open(path, "rb") as f: head = f.read(2)
data = gzip.open(path, "rb").read() if (head[0] == 31 and head[1] == 139) else open(path, "rb").read()
min_shift, depth, l_aux = struct.unpack_from("<iii", data, 4)
off = 16 + l_aux
n_ref = struct.unpack_from("<i", data, off)[0]; off += 4
print("n_ref", n_ref)
for r in range(n_ref):
    n_bin = struct.unpack_from("<i", data, off)[0]; off += 4
    print("ref", r, "n_bin", n_bin)
    for b in range(min(n_bin, 3)):
        binid, loff, n_chunk = struct.unpack_from("<IQi", data, off)
        print("  bin", binid, "loff", loff, "n_chunk", n_chunk)
        off += 16 + 16 * n_chunk
