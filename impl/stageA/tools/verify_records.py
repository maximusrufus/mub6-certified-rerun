#!/usr/bin/env python3
"""Independent cross-check of Stage A records with python-flint (arb, 160 bits).

usage: verify_records.py records.bin [N] [seed]

For N sampled records (default 12):
  1. ||F(P)||            <= froot           (arb ball upper bound)
  2. ||T^T T - I||_F     <= ort
  3. ||(I - Pi) T||_F    <= hor
  4. M-free check of the Root assertions on EXACT solutions: from P + (random horizontal
     perturbation of norm 0.04) run damped Gauss-Newton (arb arithmetic) to F ~ 0, then with
     delta = P* - P (horizontal): ||(I-TT^T) delta|| <= n0 + id s + C s^2 (+1e-9 slack for the
     finite Newton tolerance) and s^2 <= (1+ort)||T^T delta||^2 + ||(I-TT^T)delta||^2.
Exit 0 iff no failure.  Records whose Newton run does not converge are reported and skipped.
"""

import random
import struct
import sys

from flint import arb, arb_mat, ctx

ctx.prec = 160
TWO = arb(2)


def read_records(path):
    data = open(path, "rb").read()
    n = struct.unpack_from("<Q", data, 0)[0]
    off = 8
    recs = []
    for _ in range(n):
        k = struct.unpack_from("<i", data, off)[0]
        off += 4
        P = list(struct.unpack_from("<36d", data, off))
        off += 288
        T = []
        for _a in range(k):
            T.append(list(struct.unpack_from("<36d", data, off)))
            off += 288
        C, n0, idd, ort, hor, mnorm, froot, r = struct.unpack_from("<8d", data, off)
        off += 64
        recs.append(
            dict(
                k=k,
                P=P,
                T=T,
                C=C,
                n0=n0,
                id=idd,
                ort=ort,
                hor=hor,
                mnorm=mnorm,
                froot=froot,
                r=r,
            )
        )
    assert off == len(data), "trailing bytes"
    return recs


def helmert():
    H = []
    for u in range(5):
        nrm = arb((u + 1) * (u + 2)).sqrt()
        col = [
            (arb(1) / nrm if i <= u else (-arb(u + 1) / nrm if i == u + 1 else arb(0)))
            for i in range(6)
        ]
        H.append(col)
    return H  # H[u][i]


HH = helmert()


def bh_col(u, v):
    return [HH[u][i] * HH[v][j] for i in range(6) for j in range(6)]


BH = [bh_col(u, v) for u in range(5) for v in range(5)]  # 25 columns of 36


def pi_apply(x):
    rm = [sum(x[6 * i + j] for j in range(6)) / 6 for i in range(6)]
    cm = [sum(x[6 * i + j] for i in range(6)) / 6 for j in range(6)]
    tm = sum(rm) / 6
    return [x[6 * i + j] - rm[i] - cm[j] + tm for i in range(6) for j in range(6)]


def FL(P):
    F = []
    L = []
    twopi = 2 * arb.pi()
    for i in range(1, 6):
        for j in range(i):
            cs = []
            sn = []
            for g in range(6):
                ph = P[6 * i + g] - P[6 * j + g]
                cs.append((2 * ph).cos_pi())
                sn.append((2 * ph).sin_pi())
            sx = sum(cs)
            sy = sum(sn)
            F += [sy / twopi, -sx / twopi]
            rc = [arb(0) for _ in range(36)]
            rs = [arb(0) for _ in range(36)]
            for g in range(6):
                rc[6 * i + g] += cs[g]
                rc[6 * j + g] -= cs[g]
                rs[6 * i + g] += sn[g]
                rs[6 * j + g] -= sn[g]
            L += [rc, rs]
    return F, L


def norm(v):
    return sum(x * x for x in v).sqrt()


def upper(x):
    return float(x.upper())


def newton(P0, iters=60):
    P = [arb(x) for x in P0]
    for _ in range(iters):
        F, L = FL(P)
        if upper(norm(F)) < 1e-18:
            break
        # horizontal Jacobian A = L Bh (30 x 25); damped Gauss-Newton
        A = [
            [sum(L[r][p] * BH[c][p] for p in range(36)) for c in range(25)]
            for r in range(30)
        ]
        G = arb_mat(25, 25)
        for a in range(25):
            for b in range(25):
                G[a, b] = sum(A[r][a] * A[r][b] for r in range(30)) + (
                    arb("1e-9") if a == b else arb(0)
                )
        rhs = arb_mat(25, 1)
        for a in range(25):
            rhs[a, 0] = (-sum(A[r][a] * F[r] for r in range(30))).mid()
        y = G.solve(rhs)
        for p in range(36):
            P[p] = (P[p] + sum(BH[a][p] * y[a, 0] for a in range(25))).mid()
    F, _ = FL(P)
    return P, upper(norm(F))


def main():
    path = sys.argv[1]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    recs = read_records(path)
    print(f"{len(recs)} records in {path}")
    idx = list(range(len(recs)))
    rng.shuffle(idx)
    idx = idx[:N]
    bad = 0
    checked = 0
    for ii in idx:
        R = recs[ii]
        P = [arb(x) for x in R["P"]]
        F, _L = FL(P)
        fn = upper(norm(F))
        T = [[arb(x) for x in col] for col in R["T"]]
        k = R["k"]
        ortv = 0
        for a in range(k):
            for b in range(k):
                d = sum(T[a][p] * T[b][p] for p in range(36)) - (1 if a == b else 0)
                ortv += d * d
        ortv = upper(arb(ortv).sqrt()) if k else 0.0
        hv = arb(0)
        for a in range(k):
            pc = pi_apply(T[a])
            hv += sum((T[a][p] - pc[p]) ** 2 for p in range(36))
        horv = upper(hv.sqrt()) if k else 0.0
        ok = fn <= R["froot"] and ortv <= R["ort"] and horv <= R["hor"]
        msg = f"rec {ii}: k={k} |F|={fn:.2e}<=froot {R['froot']:.2e}, ort {ortv:.2e}<={R['ort']:.2e}, hor {horv:.2e}<={R['hor']:.2e}"
        # exact-solution normal check
        v = [arb(rng.gauss(0, 1)) for _ in range(25)]
        vn = norm(v)
        scale = arb(0.04) / vn
        pert = [sum(BH[a][p] * v[a] * scale for a in range(25)) for p in range(36)]
        Pn, res = newton([float((P[p] + pert[p]).mid()) for p in range(36)])
        if res < 1e-12:
            delta = [Pn[p] - P[p] for p in range(36)]
            # remove non-horizontal component caused by float start (none: perturbation horizontal)
            s = norm(delta)
            eta = [sum(T[a][p] * delta[p] for p in range(36)) for a in range(k)]
            td = [sum(T[a][p] * eta[a] for a in range(k)) for p in range(36)]
            nu = norm([delta[p] - td[p] for p in range(36)])
            sf = upper(s)
            bound = R["n0"] + R["id"] * sf + R["C"] * sf * sf + 1e-9
            c1 = upper(nu) <= bound
            etn = sum(e * e for e in eta) if k else arb(0)
            c2 = upper(s * s) <= (1 + R["ort"]) * upper(etn) + upper(nu * nu) + 1e-9
            msg += f"; exact solution at |delta|={sf:.3f}: |(I-TT^T)delta|={upper(nu):.2e} <= n_P {bound:.2e}: {c1}, square-ineq: {c2}"
            ok = ok and c1 and c2
            checked += 1
        else:
            msg += f"; Newton did not converge (res {res:.1e}), skipped"
        print(msg, "OK" if ok else "FAIL")
        bad += 0 if ok else 1
    print(
        f"summary: {len(idx)} records sampled, {checked} exact-solution normal checks, {bad} failures"
    )
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
