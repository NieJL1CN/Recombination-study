import gzip, struct, sys
path = sys.argv[1]
with open(path, "rb") as f: head = f.read(2)
if head[0] == 31 and head[1] == 139:
    data = gzip.open(path, "rb").read()
else:
    data = open(path, "rb").read()
print("size", len(data), "magic", data[:4])
min_shift, depth, l_aux = struct.unpack_from("<iii", data, 4)
print("min_shift", min_shift, "depth", depth, "l_aux", l_aux)
off = 16 + l_aux
n_ref = struct.unpack_from("<i", data, off)[0]
print("n_ref at", off, "=", n_ref)
print("aux head:", data[16:16+min(l_aux, 120)])
