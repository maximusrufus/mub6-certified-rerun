#!/usr/bin/env python3
"""Adapter check, Stage A -> Stage B hand-off (SPEC sec 4 input constraints; same bounds as stageBC validateRecord).
usage: handoff_check.py A.bin   (exit 0 iff the file parses with exact EOF and every record is within the B input bounds)"""

import struct
import sys

d = open(sys.argv[1], "rb").read()
(n,) = struct.unpack_from("<Q", d, 0)
off, bad = 8, 0
for i in range(n):
    (k,) = struct.unpack_from("<i", d, off)
    if not 0 <= k <= 10:
        print("record %d: k=%d out of range" % (i, k))
        bad += 1
        break
    P = struct.unpack_from("<36d", d, off + 4)
    T = struct.unpack_from("<%dd" % (36 * k), d, off + 4 + 288) if k else ()
    off += 4 + 288 + 288 * k
    C, n0, idd, ort, hor, mn, fr, r = struct.unpack_from("<8d", d, off)
    off += 64
    ok = (
        0 <= r <= 0.126
        and 0 <= C < 20
        and 0 <= n0 < 1e-6
        and 0 <= fr < 1e-6
        and 0 <= idd < 1e-4
        and 0 <= ort < 1e-5
        and 0 <= hor < 1e-5
        and 0 <= mn < 2
        and all(abs(x) <= 4 for x in P)
        and all(abs(x) <= 1.01 for x in T)
    )
    if not ok:
        bad += 1
        if bad < 4:
            print(
                "record %d violates B input bounds: r=%g C=%g n0=%g froot=%g id=%g ort=%g hor=%g mnorm=%g"
                % (i, r, C, n0, fr, idd, ort, hor, mn)
            )
if off != len(d):
    print("trailing/short bytes: parsed %d of %d" % (off, len(d)))
    bad += 1
print("handoff A->B: %d records, %d violations" % (n, bad))
sys.exit(1 if bad else 0)
