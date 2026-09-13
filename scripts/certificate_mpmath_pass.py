"""
30-dps certification pass (mpmath) for a sample of generators, exact coefficients parsed from the
symbolic table columns (a_p_exact = ±zeta10^k [* phi | * 1/phi | * 2], chi_p_exact = zeta5^j).
Usage: python3 certificate_mpmath_pass.py <out.tsv>
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
import csv, re, sys, math, time, numpy as np
from mpmath import mp, mpf, mpc, besselk, sqrt, sin, pi, exp, log, fabs
mp.dps = 34
N = 1951
rows = list(csv.DictReader(open(_os.path.join(_HERE,'..','data','hecke_eigenvalues_doud1951_to1e5.csv'))))
phi = (1 + sqrt(5)) / 2
z10 = [exp(2j * pi * k / 10) for k in range(10)]
z5 = [exp(2j * pi * k / 5) for k in range(5)]
def parse_a(s):
    s = s.strip()
    if s == '0': return mpc(0)
    if s == '2': return mpc(2)
    m = re.fullmatch(r'([+-]?)zeta10\^(\d+)(?:\*(.*))?', s)
    sign = -1 if m.group(1) == '-' else 1
    v = z10[int(m.group(2)) % 10] * sign
    f = m.group(3)
    if f is None: return v
    if f == 'phi': return v * phi
    if f == '1/phi': return v / phi
    if f == '2': return v * 2
    raise ValueError(s)
def parse_chi(s):
    if s == '0': return mpc(0)
    m = re.fullmatch(r'zeta5\^(\d+)', s.strip()); return z5[int(m.group(1)) % 5]
ap = {}; chi = {}
for r in rows:
    p = int(r['p']); ap[p] = parse_a(r['a_p_exact']); chi[p] = parse_chi(r['chi_p_exact'])
    # cross-check against the float columns
    assert abs(complex(ap[p]) - complex(float(r['Re_a_p']), float(r['Im_a_p']))) < 1e-12, (p, r['a_p_exact'])
primes = sorted(ap)
M = 60000
spf = np.zeros(M + 1, dtype=np.int64)
for p in primes:
    if p > M: break
    spf[p::p][spf[p::p] == 0] = p
cpk = {}
def cpow(p, k):
    if (p, k) in cpk: return cpk[(p, k)]
    if k == 0: v = mpc(1)
    elif k == 1: v = ap[p]
    elif p == N: v = ap[p] ** k
    else: v = ap[p] * cpow(p, k - 1) - chi[p] * cpow(p, k - 2)
    cpk[(p, k)] = v; return v
c = [mpc(0)] * (M + 1); c[1] = mpc(1)
for n in range(2, M + 1):
    p = int(spf[n]); m = n; k = 0
    while m % p == 0: m //= p; k += 1
    c[n] = c[m] * cpow(p, k)
# nebentypus on (Z/N)^x via primitive root 3, chi(3) = zeta5^1 (verified in the float pass)
g = 3; ind = {}; x = 1
for k in range(N - 1): ind[x] = k; x = x * g % N
def chi_of(d): return z5[ind[d % N] % 5]
def F(z, nmax):
    xx, yy = z.real, z.imag
    s = mpc(0)
    for n in range(1, nmax + 1):
        arg = 2 * pi * n * yy
        if arg > 90: break
        s += c[n] * besselk(0, arg) * sin(2 * pi * n * xx)
    return 2j * sqrt(yy) * s
def tail(y, nmax):
    return exp(-2 * pi * (nmax + 1) * y) / (1 - exp(-2 * pi * y))
gens = np.load(_os.path.join(_HERE,'..','generators','gens_unique_1951.npy'))
gens = [tuple(int(v) for v in row) for row in gens]
rng = np.random.default_rng(11)
m1 = [g for g in gens if abs(g[2]) == N]; m2 = [g for g in gens if abs(g[2]) == 2 * N]; m3 = [g for g in gens if abs(g[2]) == 3 * N]
ell = [g for g in gens if abs(g[0] + g[3]) == 1]
sample = [m1[i] for i in rng.choice(len(m1), 16, replace=False)] + [m2[i] for i in rng.choice(len(m2), 4, replace=False)] + [m3[i] for i in rng.choice(len(m3), 2, replace=False)] + ell
import os
done = set()
if os.path.exists(sys.argv[1]):
    for line in open(sys.argv[1]).read().splitlines()[1:]:
        f_ = line.split('\t'); done.add(tuple(int(v) for v in f_[:4]))
    out = open(sys.argv[1], 'a')
else:
    out = open(sys.argv[1], 'w'); out.write("a\tb\tc\td\tm\tresidual\ttail_over_scale\tsecs\n")
for (a, b, cc, d) in sample:
    if (a, b, cc, d) in done: continue
    t0 = time.time(); worst = mpf(0); fl = mpf(0)
    for s_ in (mpf('0.37'), mpf('1.23')):
        y = mpf(1) / abs(cc)
        z = mpc(mpf(-d) / cc + s_ / abs(cc), y)
        gz = (a * z + b) / (cc * z + d)
        ymin = min(y, gz.imag)
        Fz, Fg = F(z, M), F(gz, M)
        scale = max(fabs(Fz), fabs(Fg))
        R = fabs(Fg - chi_of(d) * Fz) / scale
        worst = max(worst, R); fl = max(fl, 2 * tail(ymin, M) / scale)
    line = f"{a}\t{b}\t{cc}\t{d}\t{abs(cc)//N}\t{mp.nstr(worst, 5)}\t{mp.nstr(fl, 3)}\t{time.time()-t0:.0f}\n"
    out.write(line); out.flush(); print(line, end='')
out.close()
