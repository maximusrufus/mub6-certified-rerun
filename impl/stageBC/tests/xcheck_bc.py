#!/usr/bin/env python3
"""Cross-check the Stage B / Stage C interval enclosures against python-flint (arb/acb, 220 bits).

usage: xcheck_bc.py cases.json
Every case carries exact double inputs (hex floats) and the C++ results.  For each quantity this script re-derives the
value INDEPENDENTLY from the formulas of exclusion.tex / refinement.tex with flint balls and requires:
  * Stage C quantities: the flint ball lies inside the C++ interval (exact rational comparison, as harness/flint_check.py);
  * Stage B double upper bounds: C++ value >= flint upper bound (and within 1e-6 relative slack, as a sanity check that the
    bound is not vacuous).
Exit 0 iff 0 failures.
"""

import json
import os
import sys
from fractions import Fraction

from flint import acb, arb, ctx

ctx.prec = 220
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "harness"))
try:
    from flint_check import (
        contains as _contains,
    )  # the shared exact-rational containment test
except Exception:  # pragma: no cover
    _contains = None


def hx(s):
    return float.fromhex(s)


def A(s):
    return arb(hx(s))


def frac(x):
    man, exp = x.man_exp()
    man, exp = int(man), int(exp)
    return Fraction(man * 2**exp) if exp >= 0 else Fraction(man, 2 ** (-exp))


def bounds(x):
    m, r = frac(x.mid()), frac(x.rad())
    return m - r, m + r


def contains_local(ivl, ball):
    lo, hi = Fraction(hx(ivl[0])), Fraction(hx(ivl[1]))
    bl, bu = bounds(ball)
    return lo <= bl and bu <= hi


def upper_ok(cpp, ball, slack=1e-6):
    """cpp (double) must be >= the flint upper bound, and not absurdly larger."""
    bl, bu = bounds(ball)
    v = Fraction(hx(cpp))
    return v >= bu and v <= bu * (1 + Fraction(slack)) + Fraction(1, 10**9)


def contains(ivl, ball):
    # prefer the shared exact-rational test of harness/flint_check.py when it is importable
    if _contains is not None:
        return _contains(ivl, ball)
    return contains_local(ivl, ball)


PI = arb.pi()
S6 = arb(6).sqrt()


def e(t):
    return acb(0, 2 * PI * t).exp()


def amax(lst):
    best = lst[0]
    for x in lst[1:]:
        if x.mid() > best.mid():
            best = x
    return best


def amin(a, b):
    return a if a.mid() < b.mid() else b


def rec_consts(r):
    return dict(
        k=r["k"],
        P=[A(x) for x in r["p"]],
        T=[A(x) for x in r["T"]],
        r=A(r["r"]),
        ort=A(r["ort"]),
        n0=A(r["n0"]),
        idc=A(r["id"]),
        C=A(r["C"]),
        froot=A(r["froot"]),
        mnorm=A(r["mnorm"]),
    )


def single_data(R, q, y):
    """c, J, B, S, plus the weighted error pieces for coefficient direction y (flint)."""
    S = [sum((q[j][i] for j in range(6)), acb(0)) for i in range(6)]
    a2 = [z.real**2 + z.imag**2 for z in S]
    c = [(a2[i] - 1) / (4 * PI) for i in range(6)]
    J = [[-(S[i].conjugate() * q[j][i]).imag for j in range(6)] for i in range(6)]
    T = R["T"]
    B = [
        [
            sum((J[i][j] * T[kk * 36 + 6 * j + i] for j in range(6)), arb(0))
            for kk in range(R["k"])
        ]
        for i in range(6)
    ]
    yy = y
    yinf = amax([abs(v) for v in yy])
    Sy = sum((yy[i] ** 2 * a2[i] for i in range(6)), arb(0)).sqrt()
    Syinf = amax([abs(yy[i]) * a2[i].sqrt() for i in range(6)])
    Ny = sum(
        (yy[i] ** 2 * J[i][j] ** 2 for i in range(6) for j in range(6)), arb(0)
    ).sqrt()
    Gy = amax(
        [
            abs(sum((yy[i] * S[i].conjugate() * q[j][i] for i in range(6)), acb(0)))
            for j in range(6)
        ]
    )
    return dict(S=S, a2=a2, c=c, J=J, B=B, yinf=yinf, Sy=Sy, Syinf=Syinf, Ny=Ny, Gy=Gy)


def weighted_error(R, D, X, I, use_min):
    K = arb(100001) / 100000
    nhat = R["n0"] + R["idc"] * R["r"] + R["C"] * R["r"] ** 2
    r = R["r"]
    pure = D["Gy"] * X**2
    if use_min:
        pure = amin(pure, K * D["Sy"] * I * X)
    br = (
        D["yinf"] * (K * X + r) ** 2
        + pure
        + 2 / S6 * D["Sy"] * X * r
        + 1 / S6 * D["Syinf"] * r**2
    )
    return nhat * D["Ny"] + PI * br


def check_single(c):
    R = rec_consts(c["rec"])
    M, b = c["M"], c["bin"]
    y = [A(x) for x in c["y"]]
    q = [[None] * 6 for _ in range(6)]
    for j in range(6):
        u = arb(0) if j == 0 else arb(2 * b[j - 1] + 1) / (2 * M)
        for i in range(6):
            q[j][i] = e(u - R["P"][6 * j + i]) / S6
    D = single_data(R, q, y)
    a = arb(1) / (2 * M)
    X = (arb(29) / 6).sqrt() * a
    I = arb(3) / 2 * a
    K = arb(100001) / 100000
    rH = (1 + R["ort"]).sqrt() * R["r"]
    sumg = sum(
        (
            abs(sum((y[i] * D["J"][i][j] for i in range(6)), arb(0)))
            for j in range(1, 6)
        ),
        arb(0),
    )
    bv = [sum((y[i] * D["B"][i][kk] for i in range(6)), arb(0)) for kk in range(R["k"])]
    nB = sum((v**2 for v in bv), arb(0)).sqrt()
    lhs = abs(sum((y[i] * D["c"][i] for i in range(6)), arb(0)))
    rhs = a * sumg + rH * nB + weighted_error(R, D, X, I, True)
    ok = contains(c["lhs"], lhs) and contains(c["rhs"], rhs)
    agg = sum(((D["a2"][i].sqrt() - 1) ** 2 for i in range(6)), arb(0))
    thr = (2 * PI * (K * X + R["r"])) ** 2
    tau = 2 * PI * (5 * a / S6 + R["r"])
    ok = (
        ok
        and contains(c["agg"], agg)
        and contains(c["thr"], thr)
        and contains(c["tau"], tau)
    )
    for i in range(6):
        ok = (
            ok
            and contains(c["absS2"][i], D["a2"][i])
            and contains(c["c"][i], D["c"][i])
        )
    return ok


def check_joint(c):
    R = rec_consts(c["rec"])
    r = R["r"]
    v = [[A(x) for x in c["v0"]], [A(x) for x in c["v1"]]]
    s = [A(x) for x in c["s"]]
    y = [A(x) for x in c["y"]]
    K = arb(100001) / 100000
    rH = (1 + R["ort"]).sqrt() * r
    nhat = R["n0"] + R["idc"] * r + R["C"] * r**2
    data = []
    for t in range(2):
        q = [[e(v[t][j] - R["P"][6 * j + i]) / S6 for i in range(6)] for j in range(6)]
        yt = y[6 * t : 6 * t + 6]
        data.append(single_data(R, q, yt))
    Q = [e(v[0][j] - v[1][j]) for j in range(6)]
    S = sum(Q, acb(0))
    re_ = s[0] + s[1]
    a2 = S.real**2 + S.imag**2
    absS = a2.sqrt()
    DO = [[-Q[j].imag for j in range(6)], [Q[j].real for j in range(6)]]
    DU = [(S.real * DO[0][j] + S.imag * DO[1][j]) / 6 for j in range(6)]
    nDU = sum((x**2 for x in DU), arb(0)).sqrt()
    ok = contains(c["absS"], absS)
    for isO, key in ((True, "jointO"), (False, "jointU")):
        if isO:
            ye = [y[12], y[13]]
            DE = DO
            qE = [S.real / (2 * PI), S.imag / (2 * PI)]
            EE = PI * re_**2
        else:
            ye = [y[12]]
            DE = [DU]
            qE = [(a2 / 6 - 1) / (4 * PI)]
            EE = PI * re_**2 * (1 + absS / 6)
        lhs = abs(
            sum(
                (y[i] * data[0]["c"][i] + y[6 + i] * data[1]["c"][i] for i in range(6)),
                arb(0),
            )
            + sum((ye[l] * qE[l] for l in range(len(ye))), arb(0))
        )
        rhs = arb(0)
        for t in range(2):
            sign = 1 if t == 0 else -1
            yt = y[6 * t : 6 * t + 6]
            vec = [
                sum((yt[i] * data[t]["J"][i][j] for i in range(6)), arb(0))
                + sign * sum((ye[l] * DE[l][j] for l in range(len(ye))), arb(0))
                for j in range(6)
            ]
            rhs += s[t] * sum((x**2 for x in vec), arb(0)).sqrt()
            # ball error: pure first remainder with X = s_t
            D = data[t]
            br = (
                D["yinf"] * (K * s[t] + r) ** 2
                + s[t] ** 2 * D["Gy"]
                + 2 / S6 * D["Sy"] * s[t] * r
                + 1 / S6 * D["Syinf"] * r**2
            )
            rhs += nhat * D["Ny"] + PI * br
        bv = [
            sum(
                (
                    y[i] * data[0]["B"][i][kk] + y[6 + i] * data[1]["B"][i][kk]
                    for i in range(6)
                ),
                arb(0),
            )
            for kk in range(R["k"])
        ]
        rhs += rH * sum((x**2 for x in bv), arb(0)).sqrt()
        rhs += EE * sum((x**2 for x in ye), arb(0)).sqrt()
        ok = ok and contains(c[key]["lhs"], lhs) and contains(c[key]["rhs"], rhs)
    # direct support tests
    oL = a2 / (2 * PI)
    oR = 6 * re_ * nDU + PI * re_**2 * absS
    qU = (a2 / 6 - 1) / (4 * PI)
    uR = re_ * nDU + PI * re_**2 * (1 + absS / 6)
    ok = (
        ok
        and contains(c["orthSupport"]["lhs"], oL)
        and contains(c["orthSupport"]["rhs"], oR)
    )
    ok = (
        ok
        and contains(c["unbSupport"]["lhs"], abs(qU))
        and contains(c["unbSupport"]["rhs"], uR)
    )
    return ok


def check_rel(c):
    P, Q = rec_consts(c["P"]), rec_consts(c["Q"])
    pf = [hx(x) for x in c["P"]["p"]]
    qf = [hx(x) for x in c["Q"]["p"]]

    def d0(p):
        return [
            p[6 * i + j] - p[6 * i] - p[j] + p[0] for i in range(6) for j in range(6)
        ]

    dp, dq = (
        d0(pf),
        d0(qf),
    )  # python floats (round-to-nearest) only used to pick the integer wraps
    N = [int((dq[m] - dp[m] + 0.5) // 1) for m in range(36)]

    def d0a(p):
        return [
            p[6 * i + j] - p[6 * i] - p[j] + p[0] for i in range(6) for j in range(6)
        ]

    x = [a - b - N[m] for m, (a, b) in enumerate(zip(d0a(Q["P"]), d0a(P["P"])))]
    rs = [sum((x[6 * i + j] for j in range(6)), arb(0)) for i in range(6)]
    cs = [sum((x[6 * i + j] for i in range(6)), arb(0)) for j in range(6)]
    tot = sum(rs, arb(0))
    d = [
        x[6 * i + j] - rs[i] / 6 - cs[j] / 6 + tot / 36
        for i in range(6)
        for j in range(6)
    ]
    dn = sum((t**2 for t in d), arb(0)).sqrt()
    ok = upper_ok(c["dq"], dn)
    k = P["k"]
    eta = [
        sum((P["T"][a * 36 + m] * d[m] for m in range(36)), arb(0)) for a in range(k)
    ]
    for a in range(k):
        ok = ok and contains(c["etaq"][a], eta[a])
    Ad = [
        d[m] - sum((P["T"][a * 36 + m] * eta[a] for a in range(k)), arb(0))
        for m in range(36)
    ]
    nq = sum((t**2 for t in Ad), arb(0)).sqrt()
    R_ = P["r"]
    nPR = P["n0"] + P["idc"] * R_ + P["C"] * R_**2
    D_ = R_ + dn
    v0 = nPR + nq
    lam = P["idc"] + P["C"] * D_
    e0 = P["mnorm"] * Q["froot"]
    beta = (1 + P["ort"]) * ((1 + P["ort"]) * (1 + Q["ort"])).sqrt()
    ok = (
        ok
        and upper_ok(c["D"], D_)
        and upper_ok(c["v0"], v0)
        and upper_ok(c["lam"], lam)
        and upper_ok(c["e0"], e0, 1e-3)
        and upper_ok(c["beta"], beta)
    )
    # sigma^2 >= Rayleigh quotient of G = (A S_h)^T (A S_h)  (lower bound for lambda_max)
    kq = Q["k"]
    if kq:
        W = []
        for cc in range(kq):
            col = [Q["T"][cc * 36 + m] for m in range(36)]
            rs_ = [sum((col[6 * i + j] for j in range(6)), arb(0)) for i in range(6)]
            cs_ = [sum((col[6 * i + j] for i in range(6)), arb(0)) for j in range(6)]
            t_ = sum(rs_, arb(0))
            pc = [
                col[6 * i + j] - rs_[i] / 6 - cs_[j] / 6 + t_ / 36
                for i in range(6)
                for j in range(6)
            ]
            coef = [
                sum((P["T"][a * 36 + m] * pc[m] for m in range(36)), arb(0))
                for a in range(k)
            ]
            W.append(
                [
                    pc[m]
                    - sum((P["T"][a * 36 + m] * coef[a] for a in range(k)), arb(0))
                    for m in range(36)
                ]
            )
        G = [
            [
                float(sum((W[a][m] * W[b][m] for m in range(36)), arb(0)).mid())
                for b in range(kq)
            ]
            for a in range(kq)
        ]
        vec = [1.0] * kq
        for _ in range(
            200
        ):  # power iteration (heuristic proposal of the Rayleigh vector)
            w = [sum(G[a][b] * vec[b] for b in range(kq)) for a in range(kq)]
            nrm = sum(t * t for t in w) ** 0.5
            if nrm == 0:
                break
            vec = [t / nrm for t in w]
        va = [arb(t) for t in vec]
        num = sum(
            (
                va[a] * va[b] * sum((W[a][m] * W[b][m] for m in range(36)), arb(0))
                for a in range(kq)
                for b in range(kq)
            ),
            arb(0),
        )
        den = sum((t * t for t in va), arb(0))
        lb = num / den
        sig = hx(c["sigma"])
        ok = ok and Fraction(sig) ** 2 >= bounds(lb)[0] * (1 - Fraction(1, 10**12))
    # Dhat evaluations: C++ >= exact iterate at the same number of steps
    for t, val, steps in [(hx(a), hx(b), int(hx(s))) for a, b, s in c["evals"]]:
        tt = arb(t)
        z = (1 + P["ort"]).sqrt() * tt
        b = amin(D_, (z**2 + v0**2).sqrt())
        for _ in range(steps):
            vv = amin(v0, amin((1 + P["ort"]) * b, lam * b + e0))
            psiP = (z**2 + vv**2).sqrt()
            ep = beta * tt + (arb(hx(c["sigma"]))) * vv
            nQ = Q["n0"] + Q["idc"] * b + Q["C"] * b**2
            psiQ = ((1 + Q["ort"]) * ep**2 + nQ**2).sqrt()
            b = amin(b, amin(psiP, psiQ))
        ok = ok and Fraction(val) >= bounds(b)[0] * (1 - Fraction(1, 10**12))
    return ok


def main():
    cases = json.load(open(sys.argv[1]))
    fails, counts = 0, {}
    for i, cs in enumerate(cases):
        t = cs["type"]
        counts[t] = counts.get(t, 0) + 1
        ok = {"single": check_single, "joint": check_joint, "rel": check_rel}[t](cs)
        if not ok:
            fails += 1
            if fails <= 10:
                print("CONTAINMENT FAILURE case", i, t)
    print("cases:", len(cases), counts)
    print("containment failures:", fails)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
