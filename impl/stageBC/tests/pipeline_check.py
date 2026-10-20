#!/usr/bin/env python3
"""Structural checks on a run_pipeline.sh work directory (SPEC 4.6, 4.7, 5.6).
usage: pipeline_check.py WORKDIR TAU_FINE
 * every .bin parses with exact EOF (count header, 356+288k bytes per record);
 * batches/f{j}.bin hold <=256 roots, all with written radius == TAU_FINE; fdone.txt = (#batches, #roots);
 * the concatenation of batches equals fine.bin except for the radius field;
 * every bad{j}.bin is an ORDERED SUBSEQUENCE of its batch (C exit contract).
"""

import glob
import os
import struct
import sys


def read(path):
    data = open(path, "rb").read()
    (n,) = struct.unpack_from("<Q", data, 0)
    off, recs = 8, []
    for _ in range(n):
        (k,) = struct.unpack_from("<i", data, off)
        size = 4 + 36 * 8 + k * 36 * 8 + 8 * 8
        recs.append(data[off : off + size])
        off += size
    assert off == len(data), "trailing bytes / bad size in " + path
    return recs


def radius(rec):
    return struct.unpack_from("<d", rec, len(rec) - 8)[0]


def strip_radius(rec):
    return rec[:-8]


def main():
    W, tau = sys.argv[1], float(sys.argv[2])
    fails = 0

    def check(c, msg):
        nonlocal fails
        if not c:
            fails += 1
            print("FAIL:", msg)

    fine = read(os.path.join(W, "fine.bin"))
    nb, nt = map(int, open(os.path.join(W, "batches", "fdone.txt")).read().split())
    batches = []
    for j in range(nb):
        batches.append(read(os.path.join(W, "batches", "f%d.bin" % j)))
    check(sum(len(b) for b in batches) == nt, "fdone total mismatch")
    check(nt == len(fine), "fine.bin count != fdone total")
    check(all(len(b) <= 256 for b in batches), "batch larger than 256")
    check(all(len(b) == 256 for b in batches[:-1]), "non-final batch not full")
    check(all(radius(r) == tau for b in batches for r in b), "batch radius != tau")
    cat = [strip_radius(r) for b in batches for r in b]
    check(cat == [strip_radius(r) for r in fine], "batches != fine.bin (up to radius)")
    check(all(radius(r) <= tau for r in fine), "fine.bin radius above tau")
    for j in range(nb):
        bad = read(os.path.join(W, "batches", "bad%d.bin" % j))
        it = iter(batches[j])
        ok = all(any(r == x for x in it) for r in bad)  # ordered-subsequence test
        check(ok, "bad%d.bin is not an ordered subsequence of f%d.bin" % (j, j))
    for p in glob.glob(os.path.join(W, "*.bin")):
        read(p)
    print("pipeline structure: %d batches, %d roots, %d failures" % (nb, nt, fails))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
