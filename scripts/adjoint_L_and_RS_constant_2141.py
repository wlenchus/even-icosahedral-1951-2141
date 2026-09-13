"""
The PARI half of the diagonal probe:
  (1) build L(s, Ad rho~) from the certified Galois-side table (unramified p: local poly (1-X)(1 - t_p X + X^2),
      t_p = |a_p|^2 - 2 ; p = 1951: (1 - X) ; conductor 2141^2 ; Gamma_R(s)^3 ; root number to be determined),
  (2) check its functional equation with PARI (lfuncheckfeq) for eps = +1 and eps = -1,
  (3) compute L(1, Ad) and the predicted Rankin-Selberg constant c = L(1,Ad)/zeta(2) * 2141/2142,
  (4) recompute the measured constant sum_{n<=X} |c_n|^2 / X from the table (multiplicative extension).
Everything about Ad is unconditional: L(Ad) is a Brauer ratio of Hecke L-functions (meromorphic, FE), and the
RS asymptotic follows from nonvanishing of Hecke L-functions on Re s = 1 (Wiener-Ikehara).
"""
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
import csv, math, subprocess, sys, numpy as np
from mpmath import mp, mpf, sqrt as msqrt
mp.dps = 30
N = 2141
path = _os.path.join(_HERE,'..','data','hecke_eigenvalues_doud2141_to1e5.csv')
rows = list(csv.DictReader(open(path)))
ap = {}; chi = {}; abs2 = {}
for r in rows:
    p = int(r['p'])
    ap[p] = complex(float(r['Re_a_p']), float(r['Im_a_p']))
    chi[p] = complex(float(r['Re_chi']), float(r['Im_chi']))
    a = ap[p]; abs2[p] = (a*a.conjugate()).real
primes = sorted(ap)
X = primes[-1]
print("primes:", len(primes), "max p:", X)

# ---- (4) measured RS constant from the table: c_n multiplicative, c_{p^{k+1}} = a_p c_{p^k} - chi(p) c_{p^{k-1}}; at 1951: c_{p^k} = a_N^k
def coeffs(Xmax):
    c = np.zeros(Xmax+1, dtype=complex); c[1] = 1.0
    # sieve by prime powers
    for p in primes:
        if p > Xmax: break
        pk = [1.0+0j, ap[p]]
        if p == N:
            # ramified: c_{p^k} = a_N^k
            powers = [ap[p]**k for k in range(0, 40)]
        else:
            powers = pk[:]
            while p**len(powers) <= Xmax:
                powers.append(ap[p]*powers[-1] - chi[p]*powers[-2])
        # multiply into c: process n coprime to p in increasing order via a temporary copy
        cp = c.copy()
        q = p; k = 1
        while q <= Xmax:
            # for all m with m*q <= Xmax and p ∤ m : c[m*q] += c[m] * powers[k]  (only valid if c already has all p-free parts)
            for m in range(1, Xmax//q + 1):
                if m % p: cp[m*q] += c[m]*powers[k]
            q *= p; k += 1
        c = cp
    return c
# The sieve above requires c[m] for p-free m to be complete when p is processed: process primes in increasing order,
# and since c[m] for m coprime to p only involves primes other than p (some larger!), do it differently:
def coeffs2(Xmax):
    c = np.zeros(Xmax+1, dtype=complex); c[1] = 1.0
    # compute via smallest-prime-factor factorization
    spf = np.zeros(Xmax+1, dtype=int)
    for p in primes:
        if p > Xmax: break
        for m in range(p, Xmax+1, p):
            if spf[m] == 0: spf[m] = p
    cpk = {}  # (p,k) -> c_{p^k}
    def cpow(p, k):
        if (p,k) in cpk: return cpk[(p,k)]
        if k == 0: v = 1.0+0j
        elif k == 1: v = ap[p]
        elif p == N: v = ap[p]**k
        else: v = ap[p]*cpow(p,k-1) - chi[p]*cpow(p,k-2)
        cpk[(p,k)] = v; return v
    for n in range(2, Xmax+1):
        p = spf[n]; m = n; k = 0
        while m % p == 0: m //= p; k += 1
        c[n] = c[m]*cpow(p, k)
    return c
for Xm in [10**3, 10**4, 10**5]:
    c = coeffs2(Xm)
    s2 = float(np.sum(np.abs(c[1:])**2))
    print(f"measured: sum_(n<={Xm}) |c_n|^2 / X = {s2/Xm:.5f}")
c = coeffs2(10**5)

# ---- (1)-(3) L(s, Ad) via PARI: coefficients b_n multiplicative, local poly (1-X)(1 - t X + X^2)
# b_{p^k} = h_k(1, e^{2iD}, e^{-2iD}) ; compute by power series of 1/((1-X)(1-tX+X^2))
def ad_powers(p, K):
    if p == N:
        return [1.0]*(K+1)
    t = abs2[p] - 2.0
    # 1/((1-X)(1-tX+X^2)) = sum b_k X^k ; recurrence from denominator 1 - (1+t)X + (1+t)X^2 - X^3
    b = [1.0, 1.0+t]
    b.append((1+t)*b[1] - (1+t)*b[0])
    while len(b) <= K:
        k = len(b)
        b.append((1+t)*b[k-1] - (1+t)*b[k-2] + b[k-3])
    return b[:K+1]
M = 100000
bad = np.zeros(M+1); bad[1] = 1.0
spf = np.zeros(M+1, dtype=int)
for p in primes:
    if p > M: break
    for m in range(p, M+1, p):
        if spf[m] == 0: spf[m] = p
cache = {}
for n in range(2, M+1):
    p = spf[n]; m = n; k = 0
    while m % p == 0: m //= p; k += 1
    if (p,k) not in cache:
        cache[(p,k)] = ad_powers(p, k)[k]
    bad[n] = bad[m]*cache[(p,k)]
# sanity: b_p = 1 + t_p = |a_p|^2 - 1
assert abs(bad[2] - (abs2[2]-1)) < 1e-12 and abs(bad[3] - (abs2[3]-1)) < 1e-12
_coeffs = _os.path.join(_HERE,'..','results','adjoint_coefficients_2141.gp')
with open(_coeffs, 'w') as f:
    f.write("Ad_an = [" + ",".join(f"{v:.15g}" for v in bad[1:]) + "];\n")
gp_script = r"""
default(parisize, "1G"); default(realprecision, 30);
read(COEFFS);
Ncond = 2141^2;
{
for(i=1,2,
  eps = if(i==1, 1, -1);
  L = lfuncreate([Ad_an, 0, [0,0,0], 1, Ncond, eps]);
  print("eps = ", eps, ":  lfuncheckfeq (log2 of FE residual) = ", lfuncheckfeq(L));
);
L = lfuncreate([Ad_an, 0, [0,0,0], 1, Ncond, 1]);
L1 = lfun(L, 1); print("L(1, Ad) = ", L1);
print("L(1/2, Ad) = ", lfun(L, 1/2));
print("L(2, Ad) = ", lfun(L, 2));
print("predicted RS constant L(1,Ad)/zeta(2)*2141/2142 = ", L1/zeta(2)*2141/2142);
}
"""
gp_script = gp_script.replace('read(COEFFS);', 'read("' + _coeffs + '");')
_gpfile = _os.path.join(_HERE,'..','results','adjoint_run_2141.gp')
open(_gpfile,'w').write(gp_script)
_out = _os.path.join(_HERE,'..','results','adjoint_L_gp_output_2141.txt')
subprocess.run(f'gp -q -f {_gpfile} > {_out} 2>&1', shell=True, check=False)
print(open(_out).read())
