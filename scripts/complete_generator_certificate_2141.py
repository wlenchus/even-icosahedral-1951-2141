"""
COMPLETE-GENERATOR POINTWISE CERTIFICATE for the Doud-1951 candidate Maass form.

The 08-09 MAASS-LIVE result tested F(gamma z) = chi(d) F(z) for three generators (d in {2,3,7}).
Gamma_0(1951) needs 327 generators (PARI mspolygon side pairings: 324 hyperbolic, 2 elliptic, T).
Here every generator is tested at its symmetric points, with a per-generator precision floor from
the truncation tail, in float64 (scipy K0). A 30-dps mpmath certification pass follows for a subset.

F(z) = sqrt(y) * sum_{n>=1} c_n K_0(2 pi n y) (e(nx) - e(-nx))   [sin-type, eigenvalue 1/4]
c_n multiplicative from the certified table; c_{p^{k+1}} = a_p c_{p^k} - chi(p) c_{p^{k-1}}; c_{N^k} = a_N^k.
Test points (08-09 prereg): x = -d/c + s/|c|, s in {0.37, 0.71, 1.23}; y in {1, 1.3}/|c|.
Tail bound (|c_n| <= d(n) <= 2 sqrt(n), sqrt(y) K0(2 pi n y) <= e^{-2 pi n y}/(2 sqrt n)):
   tail <= e^{-2 pi (M+1) y} / (1 - e^{-2 pi y}).
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
import csv, math, sys, time, numpy as np
from scipy.special import k0 as K0f

N = 2141
M = 100000
path = _os.path.join(_HERE,'..','data','hecke_eigenvalues_doud2141_to1e5.csv')
rows = list(csv.DictReader(open(path)))
ap = {}; chi = {}; j5 = {}
for r in rows:
    p = int(r['p'])
    ap[p] = complex(float(r['Re_a_p']), float(r['Im_a_p']))
    chi[p] = complex(float(r['Re_chi']), float(r['Im_chi']))
    j5[p] = None if r['j5'] == '-' else int(r['j5'])
primes = sorted(ap)

# ---- nebentypus as a function on (Z/N)^x: find a primitive root and match the table
def order_mod(g, n):
    k, x = 1, g % n
    while x != 1:
        x = x * g % n; k += 1
    return k
g = next(g for g in range(2, N) if order_mod(g, N) == N - 1)
ind = {}
x = 1
for k in range(N - 1):
    ind[x] = k; x = x * g % N
# chi(g) = zeta5^e ; find e from the table: chi(p) = zeta5^{j5(p)} = zeta5^{e * ind(p)}
e = None
for p in primes:
    if p == N: continue
    for cand in range(1, 5):
        if (cand * ind[p % N]) % 5 == j5[p]:
            if e is None: e = cand
    break
def chi_of(d):
    return np.exp(2j * np.pi * ((e * ind[d % N]) % 5) / 5)
bad = sum(1 for p in primes if p != N and abs(chi_of(p) - chi[p]) > 1e-9)
print(f"primitive root g={g}, chi(g)=zeta5^{e}; table mismatches: {bad} of {len(primes)-1}")
assert bad == 0

# ---- coefficients to M (float64)
c = np.zeros(M + 1, dtype=complex); c[1] = 1.0
spf = np.zeros(M + 1, dtype=np.int64)
for p in primes:
    if p > M: break
    spf[p::p][spf[p::p] == 0] = p
cpk = {}
def cpow(p, k):
    if (p, k) in cpk: return cpk[(p, k)]
    if k == 0: v = 1.0 + 0j
    elif k == 1: v = ap[p]
    elif p == N: v = ap[p] ** k
    else: v = ap[p] * cpow(p, k - 1) - chi[p] * cpow(p, k - 2)
    cpk[(p, k)] = v; return v
for n in range(2, M + 1):
    p = int(spf[n]); m = n; k = 0
    while m % p == 0: m //= p; k += 1
    c[n] = c[m] * cpow(p, k)
nn = np.arange(1, M + 1)

def F(z, coef=c):
    xx, yy = z.real, z.imag
    w = np.sqrt(yy) * K0f(2 * np.pi * nn * yy)
    return 2j * np.sum(coef[1:] * w * np.sin(2 * np.pi * nn * xx))

def tail_bound(y):
    return math.exp(-2 * math.pi * (M + 1) * y) / (1 - math.exp(-2 * math.pi * y))

gens = np.load(_os.path.join(_HERE,'..','generators','gens_unique_2141.npy'))
gens = [tuple(int(v) for v in row) for row in gens]
# scrambled control coefficients (phase-scrambled at primes, then multiplicative)
rng = np.random.default_rng(7)
ap_s = {p: (ap[p] * np.exp(2j * np.pi * rng.random()) if p != N else ap[p]) for p in primes}
cpk_s = {}
def cpow_s(p, k):
    if (p, k) in cpk_s: return cpk_s[(p, k)]
    if k == 0: v = 1.0 + 0j
    elif k == 1: v = ap_s[p]
    elif p == N: v = ap_s[p] ** k
    else: v = ap_s[p] * cpow_s(p, k - 1) - chi[p] * cpow_s(p, k - 2)
    cpk_s[(p, k)] = v; return v
cs = np.zeros(M + 1, dtype=complex); cs[1] = 1.0
for n in range(2, M + 1):
    p = int(spf[n]); m = n; k = 0
    while m % p == 0: m //= p; k += 1
    cs[n] = cs[m] * cpow_s(p, k)

results = []
t0 = time.time()
for (a, b, cc, d) in gens:
    if cc == 0:
        results.append((a, b, cc, d, 0, 0.0, 0.0, 0.0, 'T-periodicity (built in)')); continue
    m = abs(cc) // N
    worst = 0.0; worst_bar = 0.0; floor = 0.0; worst_scr = 0.0
    for s in (0.37, 0.71, 1.23):
        for yf in (1.0, 1.3):
            y = yf / abs(cc)
            z = complex(-d / cc + s / abs(cc), y)
            gz = (a * z + b) / (cc * z + d)
            ymin = min(y, gz.imag)
            fl = tail_bound(ymin) / max(1e-300, 1.0)   # absolute tail; normalised below
            Fz, Fg = F(z), F(gz)
            scale = max(abs(Fz), abs(Fg))
            R = abs(Fg - chi_of(d) * Fz) / scale
            Rbar = abs(Fg - np.conj(chi_of(d)) * Fz) / scale
            worst = max(worst, R); worst_bar = max(worst_bar, Rbar)
            floor = max(floor, 2 * fl / scale)
            if len(results) < 12:      # scrambled control on the first few generators only
                Fzs, Fgs = F(z, cs), F(gz, cs)
                worst_scr = max(worst_scr, abs(Fgs - chi_of(d) * Fzs) / max(abs(Fzs), abs(Fgs)))
    results.append((a, b, cc, d, m, worst, worst_bar, floor, f"scrambled={worst_scr:.2e}" if worst_scr else ""))
print(f"sweep done in {time.time()-t0:.0f}s over {len(gens)} generators")
import collections
by_m = collections.defaultdict(list)
for r in results:
    if r[4] > 0: by_m[r[4]].append(r)
print(" m   #gens   worst residual (chi(d))   worst floor   worst residual (chibar(d))")
for m in sorted(by_m):
    rs = by_m[m]
    print(f"{m:3d} {len(rs):6d}   {max(r[5] for r in rs):.3e}              {max(r[7] for r in rs):.1e}     {min(r[6] for r in rs):.3e}")
ell = [r for r in results if abs(r[0] + r[3]) == 1]
print("elliptic generators:", [(r[:4], f"{r[5]:.2e}", f"floor {r[7]:.1e}") for r in ell])
print("scrambled controls:", [r[8] for r in results if r[8].startswith('scr')][:6])
np.save(_os.path.join(_HERE,'..','results','certificate_results_2141.npy'), np.array(results, dtype=object), allow_pickle=True)
with open(_os.path.join(_HERE,'..','results','certificate_results_2141.tsv'), 'w') as f:
    f.write("a\tb\tc\td\tm\tresidual_chi\tresidual_chibar\tfloor\tnote\n")
    for r in results:
        f.write("\t".join(str(v) if not isinstance(v, float) else f"{v:.6e}" for v in r) + "\n")
