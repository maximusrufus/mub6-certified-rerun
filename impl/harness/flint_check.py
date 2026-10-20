#!/usr/bin/env python3
"""Re-evaluate C++ interval records with python-flint acb (prec >= 128) and check CONTAINMENT.

usage: flint_check.py records.json [--prec 160]
Input: JSON list written by tests/gen_records.cpp (hex-float strings => exact doubles).
Each record carries the C++ enclosure; flint evaluates the same quantity at the recorded sample
points (inside the input intervals) and the C++ enclosure must contain the flint ball, compared
EXACTLY with rational arithmetic. Exit 0 iff 0 failures.
"""

import json
import sys
from fractions import Fraction

from flint import acb, arb, ctx

PREC = 160


def hx(s):
    return float.fromhex(s)


def _frac(x):
    """exact Fraction of an arb with exact midpoint (ignores its radius)."""
    man, exp = x.man_exp()
    man, exp = int(man), int(exp)
    return Fraction(man * 2 ** exp) if exp >= 0 else Fraction(man, 2 ** (-exp))


def ball_bounds(x):
    """exact rational [lower, upper] of an arb ball."""
    fm, fr = _frac(x.mid()), _frac(x.rad())
    return fm - fr, fm + fr


def contains(iv_hex, x):
    lo, hi = Fraction(hx(iv_hex[0])), Fraction(hx(iv_hex[1]))
    bl, bu = ball_bounds(x)
    if lo == hi:
        # Degenerate (exactly representable) C++ result: an exact claim such as cos(pi/2) = 0.
        # A flint ball of radius ~2^-prec around the exact value cannot be contained in a point, so the
        # rigorous test is that the point lies INSIDE the ball (a disjoint ball refutes the claim).
        return bl <= lo <= bu
    return lo <= bl and bu <= hi


def A(x):
    return arb(hx(x))


def cont_c(civ, z):
    return contains(civ["re"], z.real) and contains(civ["im"], z.imag)


def cz(c):
    return acb(A(c["re"][2]), A(c["im"][2]))


def cvec(v):
    return [cz(c) for c in v]


def check(rec):
    op = rec["op"]
    if op in ("add", "sub", "mul", "div"):
        a, b = A(rec["a"][2]), A(rec["b"][2])
        r = {"add": a + b, "sub": a - b, "mul": a * b, "div": a / b}[op]
        return contains(rec["out"], r)
    if op == "sqr":
        a = A(rec["a"][2])
        return contains(rec["out"], a * a)
    if op == "sqrt":
        return contains(rec["out"], A(rec["a"][2]).sqrt())
    if op in ("cmul", "cdiv"):
        a, b = cz(rec["a"]), cz(rec["b"])
        return cont_c(rec["out"], a * b if op == "cmul" else a / b)
    if op == "cabs2":
        a = cz(rec["a"])
        return contains(rec["out"], a.real * a.real + a.imag * a.imag)
    if op == "ip":
        u, v = cvec(rec["u"]), cvec(rec["v"])
        ip = sum((x.conjugate() * y for x, y in zip(u, v)), acb(0))
        a2 = ip.real * ip.real + ip.imag * ip.imag
        n2 = sum((x.real * x.real + x.imag * x.imag for x in u), arb(0))
        return (
            cont_c(rec["ip"], ip)
            and contains(rec["abs2"], a2)
            and contains(rec["resid"], a2 - arb(1) / 6)
            and contains(rec["norm2u"], n2)
        )
    if op == "unit_root":
        z = (acb(0, 2) * arb.pi() * arb(rec["k"]) / rec["N"]).exp()
        return cont_c(rec["out"], z)
    if op == "consts":
        return (
            contains(rec["one_sixth"], arb(1) / 6)
            and contains(rec["inv_sqrt6"], 1 / arb(6).sqrt())
            and contains(rec["sqrt6"], arb(6).sqrt())
            and contains(rec["sqrt3"], arb(3).sqrt())
            and contains(rec["pi"], arb.pi())
        )
    if op == "mub":
        return check_mub(rec)
    raise ValueError("unknown op " + op)


def flint_bases():
    I = acb(0, 1)
    s6 = arb(6).sqrt()
    std = [[acb(1 if j == k else 0) for j in range(6)] for k in range(6)]
    w = lambda m: (I * 2 * arb.pi() * m / 6).exp()
    four = [[w(j * k) / s6 for j in range(6)] for k in range(6)]
    s3 = arb(3).sqrt()
    d = acb((1 - s3) / 2, (s3 / 2).sqrt())
    dc = d.conjugate()
    c = [acb(1), I * d, -d, -I, -dc, I * dc]
    bj = [[c[(k - j) % 6] / s6 for j in range(6)] for k in range(6)]
    return [std, four, bj]


def check_mub(rec):
    B = flint_bases()
    ok = True
    for b in range(3):
        for k in range(6):
            for j in range(6):
                if not cont_c(rec["bases"][b][k][j], B[b][k][j]):
                    ok = False
    pairs = [(0, 1), (0, 2), (1, 2)]
    for p, (x, y) in enumerate(pairs):
        for a in range(6):
            for b in range(6):
                ip = sum(
                    (B[x][a][j].conjugate() * B[y][b][j] for j in range(6)), acb(0)
                )
                a2 = ip.real * ip.real + ip.imag * ip.imag
                if not contains(rec["abs2"][p][6 * a + b], a2):
                    ok = False
                # Flint truth check on the literature triple: |<u,v>|^2 = 1/6 to prec (sanity of the construction)
                if abs(float((a2 - arb(1) / 6).mid())) > 1e-30:
                    print("flint: pair", (x, y), a, b, "not unbiased:", a2)
                    ok = False
    return ok


def main():
    global PREC
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if "--prec" in sys.argv:
        PREC = int(sys.argv[sys.argv.index("--prec") + 1])
    ctx.prec = max(PREC, 128)
    recs = json.load(open(args[0]))
    fails, counts = 0, {}
    for i, r in enumerate(recs):
        counts[r["op"]] = counts.get(r["op"], 0) + 1
        if not check(r):
            fails += 1
            if fails <= 10:
                print("CONTAINMENT FAILURE record", i, json.dumps(r)[:400])
    print("records:", len(recs), counts)
    print("containment failures:", fails)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
