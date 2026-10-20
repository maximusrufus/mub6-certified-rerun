#!/usr/bin/env python3
"""S6: independent python-flint (arb/acb, 200 bit) re-verification of a seeded random sample of EXCLUDED Stage A cells.

usage: verify_excluded.py seeds.txt excl_sample.txt [procs]

excl_sample.txt is written by `arun --sample-excl M file` (a deterministic splitmix64 hash sample, 1/M of all excluded
cells of both kinds).  For each sampled cell the seed context (frames, D-span, variation matrices), the cell matrix D,
and the exclusion inequality are RE-DERIVED here from the seed phases q only (nothing but the cell index, the claimed
constants h, z2, Er, Rs and, for directional exclusions, the exact certificate direction y[20] is taken from the C++ run):
  kind 0 (norm):   lower( sum_ij (|D_ij| - 1)^2 ) > z2_claimed
  kind 1 (direct): lower( y.c ) > upper( Rs|g_seed| + h sum|g_ang| + Er sqrt(sum D2 y^2) + z2/2 max(-y_D)+ + Rs^2/2 (|yX|+|yY|) )
and each claimed constant must be >= its independently recomputed rigorous upper bound (Rs, h, z2, Er).
A "disagreement" is any sampled cell for which an inequality or a claim check fails (or the re-derivation fails).
Prints "disagreements: N" and exits 0 iff N == 0 and at least one cell was checked.
"""

import math
import sys
from multiprocessing import Pool

from flint import acb, acb_mat, arb, ctx

ctx.prec = 200
RS_CYCLES = 0.0320001  # Params::rsCycles
PI = arb.pi()
I_ = acb(0, 1)


def cis_turns(t):  # e(t) = cos 2 pi t + i sin 2 pi t
    x = 2 * t
    return acb(x.cos_pi(), x.sin_pi())


def fl(
    x,
):  # arb -> python float midpoint (only used for the C++-mirrored heuristic argmax choices)
    return float(x.mid().str(20, radius=False))


def mat(rows):
    return acb_mat([[acb(z) for z in r] for r in rows])


def adj(A):
    return A.transpose().conjugate()


def hull(a, b):
    return a.union(b)


def frame_k(p, bk, bl):
    S = p[0] + p[1] + p[2] + p[3]
    Sabs = abs(S)
    half = acb(arb(1) / 2)
    q1 = [half] * 4
    c1 = S / 2
    r2 = [p[k] - q1[k] * c1 for k in range(4)]
    n2 = sum((abs(z) ** 2 for z in r2), arb(0)).sqrt()
    if not n2 > arb("0.449999"):
        raise ValueError("frame: GS residual 2 small")
    q2 = [z / n2 for z in r2]
    res, nrm = [], []
    for k in range(4):
        v = []
        for j in range(4):
            e = acb(1 if j == k else 0)
            v.append(e - q1[j] * half - q2[j] * q2[k].conjugate())
        res.append(v)
        nrm.append(sum((abs(z) ** 2 for z in v), arb(0)).sqrt())
    if not nrm[bk] > arb("0.449999"):
        raise ValueError("frame: GS residual 3 small")
    q3 = [z / nrm[bk] for z in res[bk]]
    res2, nrm2 = {}, {}
    for l in range(4):
        if l == bk:
            continue
        v = []
        for j in range(4):
            e = acb(1 if j == l else 0)
            v.append(
                e - q1[j] * half - q2[j] * q2[l].conjugate() - q3[j] * q3[l].conjugate()
            )
        res2[l] = v
        nrm2[l] = sum((abs(z) ** 2 for z in v), arb(0)).sqrt()
    if not nrm2[bl] > arb("0.449999"):
        raise ValueError("frame: GS residual 4 small")
    q4 = [z / nrm2[bl] for z in res2[bl]]
    K = mat([[q3[j], q4[j]] for j in range(4)])
    return K, Sabs


class Seed:
    def __init__(self, q, fc):
        self.q = q
        qa = [arb(x) for x in q]  # exact doubles
        self.a = cis_turns(qa[0])
        self.b = [cis_turns(qa[1 + j]) for j in range(4)]
        self.c = [cis_turns(qa[5 + i]) for i in range(4)]
        one = acb(1)
        self.Rs_true = 2 * PI * arb(math.nextafter(RS_CYCLES, 1.0))
        E0 = mat([[one, one], [one, self.a]])
        B0 = mat([[one] * 4, self.b])
        C0 = mat([[one, self.c[i]] for i in range(4)])
        KB, SbAbs = frame_k([z.conjugate() for z in self.b], fc[0], fc[1])
        KC, ScAbs = frame_k(list(self.c), fc[2], fc[3])
        self.mB = (arb(4) - SbAbs).sqrt()
        self.mC = (arb(4) - ScAbs).sqrt()
        # imin: only the lower end of m is used downstream (m divides); use it as an exact point (no hull inflation)
        self.m = self.mB.lower() if self.mB.lower() < self.mC.lower() else self.mC.lower()
        if not self.m > arb("1.24"):
            raise ValueError("frame: m too small")
        self.Cop = (arb(4) + ScAbs).sqrt()
        B0s = adj(B0)
        GB = B0s * (B0 * B0s).inv()
        C0s = adj(C0)
        GC = C0 * (C0s * C0).inv()
        KBs = adj(KB)
        s6 = arb(6).sqrt()
        self.s6 = s6
        self.KC6 = mat([[KC[i, j] * s6 for j in range(2)] for i in range(4)])
        self.KBs = KBs
        self.PK = KB * KBs
        GBs = adj(GB)
        E0s = adj(E0)
        self.Dspan = -(C0 * E0s * GBs)
        self.A0 = (arb(2) + abs(one + self.a)).sqrt()
        self.Cst, self.Lft, self.Lc = [], [], []
        for l in range(9):
            E1 = acb_mat(2, 2)
            B1 = acb_mat(2, 4)
            C1 = acb_mat(4, 2)
            if l == 0:
                E1[1, 1] = I_ * self.a
            elif l <= 4:
                B1[1, l - 1] = I_ * self.b[l - 1]
            else:
                C1[l - 5, 1] = I_ * self.c[l - 5]
            E1s, B1s, C1s = adj(E1), adj(B1), adj(C1)
            t1 = C1 * E0s * GBs
            t2 = C0 * E1s * GBs
            t3 = GC * E0s * B1 * self.PK
            self.Cst.append(-(t1 + t2 + t3))
            self.Lft.append(B1s * GBs)
            self.Lc.append(GC * C1s)
        sx = one + self.a + sum(self.b, acb(0))
        sy = one + self.a + sum(self.c, acb(0))
        self.cX = [sx.imag, -sx.real]
        self.cY = [sy.imag, -sy.real]

    def level_bounds(self, Nu, Rs_claim, z2_claim):
        """rigorous (ball) h, p, z2, Er recomputed with the claimed Rs (point) and, for Er, the claimed z2 (point)."""
        h = PI / Nu
        Rs = arb(Rs_claim)
        t = self.A0 + Rs
        Astar = arb.union(min_ball(t.lower(), arb(2)), min_ball(t.upper(), arb(2)))
        qq = Rs**2 * (1 + (Astar**2 + 6) / self.mC**2)
        p = qq / (self.s6 + (arb(6) - qq).sqrt())
        off = Rs * (self.s6 + Astar) / self.m
        z2 = off**2 + (6 * h + p) ** 2
        zz = arb(z2_claim).sqrt()
        Lx = self.Cop * (1 + (self.A0 / self.mB) ** 2).sqrt()
        Lbig = (Lx.upper() if Lx.upper() > self.s6.upper() else self.s6.upper())  # imax: only the upper end is used
        eoff = Rs**2 * (Lbig + self.A0 + 1) / (self.m * 2) + zz * Rs / self.m
        e3 = arb(252).sqrt() / 2 * h**2 + p
        Er = (eoff**2 + e3**2).sqrt()
        return h, z2, Er


def min_ball(x, y):
    return min(y, x)


def cell_D(S, Nu, ia, it, ib, ig):
    def mid(i):
        return arb(2 * i + 1) / (2 * Nu)

    eA, eB, eG, eT = (
        cis_turns(mid(ia)),
        cis_turns(mid(ib)),
        cis_turns(mid(ig)),
        cis_turns(mid(it)),
    )
    cT, sT = eT.real, eT.imag
    ab, ag = eA * eB, eA * eG
    abc, agc = eA * eB.conjugate(), eA * eG.conjugate()
    U = mat([[ab * cT, ag * sT], [-(agc * sT), abc * cT]])
    D = S.Dspan + S.KC6 * U * S.KBs
    return U, D, (eA, eB, eG, cT, sT, ab, ag, abc, agc)


def check_cell(S, rec):
    kind, Nu, ia, it, ib, ig, h_c, z2_c, Er_c, Rs_c, y, fc = rec
    # claims >= independently recomputed rigorous upper bounds
    if not arb(Rs_c) >= S.Rs_true:
        return "Rs claim < 2 pi rs"
    h, z2, Er = S.level_bounds(Nu, Rs_c, z2_c)
    if not arb(h_c) >= h:
        return "h claim < pi/Nu"
    if not arb(z2_c) >= z2:
        return "z2 claim < recomputed z2"
    if not arb(Er_c) >= Er:
        return "Er claim < recomputed Er"
    U, D, (eA, eB, eG, cT, sT, ab, ag, abc, agc) = cell_D(S, Nu, ia, it, ib, ig)
    D2 = [abs(D[i, j]) ** 2 for i in range(4) for j in range(4)]
    if kind == 0:
        s = sum(((abs(D[i, j]) - 1) ** 2 for i in range(4) for j in range(4)), arb(0))
        if not s.lower() > arb(z2_c):
            return "norm exclusion not reproduced (lower=%s z2=%s)" % (
                s.lower().str(8),
                z2_c,
            )
        return None
    # directional: Jacobian columns and the certificate inequality
    DPK = D * S.PK
    Jc = []
    for l in range(9):
        t = S.Cst[l] - D * S.Lft[l]
        Jc.append(t - S.Lc[l] * DPK)
    dU = [acb_mat(2, 2) for _ in range(4)]
    for i in range(2):
        for j in range(2):
            dU[0][i, j] = I_ * U[i, j]
    dU[1][0, 0] = -(ab * sT)
    dU[1][0, 1] = ag * cT
    dU[1][1, 0] = -(agc * cT)
    dU[1][1, 1] = -(abc * sT)
    dU[2][0, 0] = I_ * U[0, 0]
    dU[2][1, 1] = -(I_ * U[1, 1])
    dU[3][0, 1] = I_ * U[0, 1]
    dU[3][1, 0] = -(I_ * U[1, 0])
    for a in range(4):
        Jc.append(S.KC6 * dU[a] * S.KBs)
    J = [[arb(0)] * 13 for _ in range(20)]
    J[0][0], J[1][0], J[2][0], J[3][0] = S.a.real, S.a.imag, S.a.real, S.a.imag
    for j in range(4):
        J[0][1 + j], J[1][1 + j] = S.b[j].real, S.b[j].imag
        J[2][5 + j], J[3][5 + j] = S.c[j].real, S.c[j].imag
    c = [S.cX[0], S.cX[1], S.cY[0], S.cY[1]] + [arb(0)] * 16
    for i in range(4):
        for j in range(4):
            f = 4 * i + j
            c[4 + f] = (D2[f] - 1) / 2
            for l in range(13):
                Jl = Jc[l][i, j]
                J[4 + f][l] = D[i, j].real * Jl.real + D[i, j].imag * Jl.imag
    yv = [arb(v) for v in y]
    lhs = sum((yv[i] * c[i] for i in range(20)), arb(0))
    g = [sum((yv[i] * J[i][l] for i in range(20)), arb(0)) for l in range(13)]
    Rs, hh, z2c, Erc = arb(Rs_c), arb(h_c), arb(z2_c), arb(Er_c)
    rhs = Rs * sum((g[l] ** 2 for l in range(9)), arb(0)).sqrt()
    rhs += hh * sum((abs(g[l]) for l in range(9, 13)), arb(0))
    rhs += Erc * sum((D2[f] * yv[4 + f] ** 2 for f in range(16)), arb(0)).sqrt()
    mx = arb(0)
    for f in range(16):
        mx = max(mx, -yv[4 + f])
    rhs += z2c / 2 * mx
    yx = (yv[0] ** 2 + yv[1] ** 2).sqrt()
    yy = (yv[2] ** 2 + yv[3] ** 2).sqrt()
    rhs += Rs**2 / 2 * (yx + yy)
    if not lhs.lower() > rhs.upper():
        return "directional exclusion not reproduced (lhs.lo=%s rhs.hi=%s)" % (
            lhs.lower().str(8),
            rhs.upper().str(8),
        )
    return None


def work(args):
    idx, qd, recs = args
    fc = recs[0][11]  # frame choices (same for all cells of a seed)
    out = []
    try:
        S = Seed(qd, fc)
    except Exception as e:  # noqa
        return [(idx, len(recs), "seed re-derivation failed: %s" % e)]
    for r in recs:
        try:
            msg = check_cell(S, r)
        except Exception as e:  # noqa
            msg = "exception: %r" % (e,)
        out.append((idx, 1, msg))
    return out


def main():
    seeds_path, ex_path = sys.argv[1], sys.argv[2]
    procs = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    seeds = []
    for ln in open(seeds_path):
        if not ln.strip() or ln.startswith("#"):
            continue
        t = ln.split()
        seeds.append(
            [((float(b) + 0.5) / 112.0) for b in t[1:10]]
        )  # same binary64 division as arun::parseLine
    by = {}
    n = 0
    for ln in open(ex_path):
        if ln.startswith("#") or not ln.strip():
            continue
        t = ln.split()
        idx, kind, Nu, ia, it, ib, ig = (int(x) for x in t[:7])
        h, z2, Er, Rs = (float.fromhex(x) for x in t[7:11])
        fc = [int(ch) for ch in t[11]]
        y = [float.fromhex(x) for x in t[12:32]] if kind == 1 else None
        by.setdefault(idx, []).append((kind, Nu, ia, it, ib, ig, h, z2, Er, Rs, y, fc))
        n += 1
    jobs = [(i, seeds[i], r) for i, r in sorted(by.items())]
    bad = []
    checked = {0: 0, 1: 0}
    with Pool(procs) as pool:
        for res in pool.imap_unordered(work, jobs, chunksize=4):
            for idx, cnt, msg in res:
                if msg:
                    bad.append((idx, msg))
    for idx, r in by.items():
        for x in r:
            checked[x[0]] += 1
    for idx, msg in bad[:20]:
        print("DISAGREE seed %d: %s" % (idx, msg))
    print(
        "sampled cells checked: %d (norm %d, directional %d) over %d seeds"
        % (n, checked[0], checked[1], len(by))
    )
    print("disagreements: %d" % len(bad))
    sys.exit(0 if (not bad and n > 0) else 1)


if __name__ == "__main__":
    main()
