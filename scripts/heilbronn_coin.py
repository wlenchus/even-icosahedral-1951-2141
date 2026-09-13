"""
Heilbronn coin for the binary icosahedral group 2.A5 (the even-icosahedral core).

At a point s0 != 1, the Heilbronn character theta = sum_chi n_chi chi with n_chi = ord_{s0} L(s, chi).
Known positivity (all L-functions below are entire away from s = 1):
  (H) Hecke:      <theta, Ind_H lambda> >= 0  for every subgroup H and 1-dim lambda of H
  (L) Langlands:  <theta, Ind_{2.A4} (2|_{2.A4} (x) omega^j)> >= 0  (tetrahedral base change; all proper
                  subgroups of 2.A5 are solvable, and the only non-monomial 2-dim irreps of proper subgroups
                  are the three faithful 2-dim irreps of 2.A4)
  (Z) H = 1:      N := theta(1) = ord_{s0} zeta_L >= 0
Artin for chi at s0  <=>  n_chi >= 0.  Question: the minimal N compatible with n_2 = -1 (a simple pole of
L(rho~) at s0), and the configurations achieving it.  Everything is exact integer arithmetic on class
functions computed numerically from unit icosians (verified by orthogonality).
"""
import itertools, math, sys
import numpy as np

phi = (1 + 5 ** 0.5) / 2
# ---------- 1. the 120 unit icosians ----------
els = []
for i in range(4):
    for s in (1, -1):
        v = [0.0] * 4; v[i] = s; els.append(tuple(v))
for signs in itertools.product((0.5, -0.5), repeat=4):
    els.append(tuple(signs))
even_perms = [p for p in itertools.permutations(range(4))
              if sum(1 for a in range(4) for b in range(a + 1, 4) if p[a] > p[b]) % 2 == 0]
base = (0.0, 0.5, 0.5 / phi, 0.5 * phi)
for p in even_perms:
    for signs in itertools.product((1, -1), repeat=3):
        vals = [0.0, signs[0] * base[1], signs[1] * base[2], signs[2] * base[3]]
        v = [0.0] * 4
        for k in range(4):
            v[p[k]] = vals[k]
        els.append(tuple(v))
els = sorted(set(tuple(round(c, 12) for c in v) for v in els))
assert len(els) == 120, len(els)
E = np.array(els)

def qmul(a, b):
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    return (a0*b0 - a1*b1 - a2*b2 - a3*b3,
            a0*b1 + a1*b0 + a2*b3 - a3*b2,
            a0*b2 - a1*b3 + a2*b0 + a3*b1,
            a0*b3 + a1*b2 - a2*b1 + a3*b0)

def find(v):
    d = np.abs(E - np.array(v)).sum(axis=1)
    k = int(np.argmin(d)); assert d[k] < 1e-6, d[k]
    return k

n = 120
M = np.zeros((n, n), dtype=int)
for i in range(n):
    for j in range(n):
        M[i, j] = find(qmul(els[i], els[j]))
e = find((1, 0, 0, 0))
inv = np.zeros(n, dtype=int)
for i in range(n):
    inv[i] = int(np.where(M[i] == e)[0][0])
# orders
def order(g):
    k, x = 1, g
    while x != e:
        x = M[x, g]; k += 1
    return k
ords = np.array([order(g) for g in range(n)])
# conjugacy classes
cls = -np.ones(n, dtype=int); reps = []
for g in range(n):
    if cls[g] >= 0: continue
    c = len(reps); reps.append(g)
    for x in range(n):
        cls[M[M[x, g], inv[x]]] = c
ncls = len(reps)
csize = np.array([(cls == c).sum() for c in range(ncls)])
print("classes:", ncls, "sizes:", csize.tolist(), "orders:", [int(ords[r]) for r in reps])

# ---------- 2. irreducible characters ----------
q0 = E[:, 0]
def conj_gal(x):
    # sqrt5 -> -sqrt5 on the possible values of q0
    table = {0.0: 0.0, 1.0: 1.0, -1.0: -1.0, 0.5: 0.5, -0.5: -0.5,
             round(phi/2, 12): round(-1/(2*phi), 12), round(-phi/2, 12): round(1/(2*phi), 12),
             round(1/(2*phi), 12): round(-phi/2, 12), round(-1/(2*phi), 12): round(phi/2, 12)}
    return table[round(float(x), 12)]
q0p = np.array([conj_gal(x) for x in q0])
def U(k, x):
    # Chebyshev U_k(x) = sin((k+1)t)/sin t, x = cos t ; character of Sym^k of SU(2)
    u0, u1 = np.ones_like(x), 2 * x
    if k == 0: return u0
    for _ in range(k - 1):
        u0, u1 = u1, 2 * x * u1 - u0
    return u1
chars = {
    '1':  np.ones(n),
    '2':  2 * q0,
    "2'": 2 * q0p,
    '3':  U(2, q0),
    "3'": U(2, q0p),
    '4':  (2 * q0) * (2 * q0p),
    "4'": U(3, q0),
    '5':  U(4, q0),
    '6':  U(5, q0),
}
names = list(chars.keys())
dims = {k: int(round(chars[k][e])) for k in names}
def ip(a, b):
    return float(np.dot(a, b)) / n   # real characters here
G = np.array([[ip(chars[a], chars[b]) for b in names] for a in names])
assert np.allclose(G, np.eye(len(names)), atol=1e-9), G
assert np.allclose(U(3, q0p), chars["4'"]) and np.allclose(U(5, q0p), chars['6'])
print("irreps orthonormal; dims:", dims)
# class-function values table for the record
print("character table (classes ordered as above):")
for k in names:
    print(f"  {k:>3}: ", [round(float(chars[k][r]), 4) for r in reps])

# ---------- 3. subgroups up to conjugacy ----------
def closure(gens):
    S = {e}; frontier = [e]
    gens = list(gens)
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = M[x, g]
                if y not in S:
                    S.add(y); new.append(y)
        frontier = new
    return frozenset(S)
subs = set()
for a in range(n):
    for b in range(a, n):
        subs.add(closure((a, b)))
subs = list(subs)
print("subgroups (all):", len(subs))
def conjugate_set(H, x):
    return frozenset(M[M[x, h], inv[x]] for h in H)
classes_of_subs = []
seen = set()
for H in sorted(subs, key=len):
    if H in seen: continue
    orbit = {conjugate_set(H, x) for x in range(n)}
    seen |= orbit
    classes_of_subs.append(H)
print("subgroup classes:", [(len(H)) for H in classes_of_subs])

# ---------- 4. one-dimensional characters of a subgroup (brute force homomorphisms) ----------
def one_dim_chars(H):
    H = sorted(H)
    m = 1
    for h in H:
        m = m * ords[h] // math.gcd(m, ords[h])   # exponent
    # generating pairs
    gens = None
    for a in H:
        if closure((a,)) == frozenset(H): gens = (a,); break
    if gens is None:
        for a in H:
            for b in H:
                if closure((a, b)) == frozenset(H): gens = (a, b); break
            if gens: break
    assert gens is not None, "need >2 generators"
    out = []
    for vals in itertools.product(range(m), repeat=len(gens)):
        lam = {e: 0}; ok = True
        frontier = [e]
        while frontier and ok:
            new = []
            for x in frontier:
                for g, v in zip(gens, vals):
                    y = M[x, g]; val = (lam[x] + v) % m
                    if y in lam:
                        if lam[y] != val: ok = False; break
                    else:
                        lam[y] = val; new.append(y)
                if not ok: break
            frontier = new
        if ok and len(lam) == len(H):
            out.append({h: np.exp(2j * np.pi * lam[h] / m) for h in H})
    # dedupe
    uniq = []
    for lam in out:
        key = tuple(np.round([lam[h] for h in H], 9))
        if key not in {tuple(np.round([u[h] for h in H], 9)) for u in uniq}:
            uniq.append(lam)
    return uniq

def induce(H, tau):
    """tau: dict h -> complex value (a character of H). Returns class function on G (complex)."""
    Hs = set(H)
    out = np.zeros(n, dtype=complex)
    for g in range(n):
        s = 0
        for x in range(n):
            y = M[M[inv[x], g], x]
            if y in Hs: s += tau[y]
        out[g] = s / len(H)
    return out

def decompose(cf):
    return np.array([np.dot(chars[k], np.conj(cf)).real / n for k in names])

rows = []   # (label, multiplicity vector)
for H in classes_of_subs:
    for lam in one_dim_chars(H):
        v = decompose(induce(H, lam))
        assert np.allclose(v, np.round(v), atol=1e-8) and (np.round(v) >= 0).all(), v
        rows.append((f"Ind_H{len(H)} lambda", np.round(v).astype(int)))
# tetrahedral: 2.A4 is the order-24 subgroup class
H24 = [H for H in classes_of_subs if len(H) == 24][0]
for lam in one_dim_chars(H24):
    tau = {h: 2 * q0[h] * lam[h] for h in H24}      # 2|_{2.A4} (x) omega^j
    v = decompose(induce(H24, tau))
    assert np.allclose(v, np.round(v), atol=1e-8), v
    rows.append(("Ind_{2.A4}(2 (x) omega^j)", np.round(v).astype(int)))
# dedupe rows
uniq = {}
for lab, v in rows:
    uniq.setdefault(tuple(v), lab)
rows = [(lab, np.array(v)) for v, lab in uniq.items()]
print("distinct positivity constraints:", len(rows))
for lab, v in rows:
    print("  ", lab, dict(zip(names, v.tolist())))

# ---------- 5. minimal pole configurations ----------
A = np.array([v for _, v in rows])          # constraints A n >= 0
d = np.array([dims[k] for k in names])
def minimal_configs(target, Nmax=8):
    ti = names.index(target)
    best = None; sols = []
    for N in range(1, Nmax + 1):
        # brute force: n_target = -1, others in [-N, N], sum d n = N, A n >= 0
        free = [i for i in range(len(names)) if i != ti]
        rng = range(-N, N + 1)
        found = []
        # prune with the H=1 and Z rows via itertools is fine at N<=4
        for vals in itertools.product(rng, repeat=len(free)):
            nvec = np.zeros(len(names), dtype=int); nvec[ti] = -1
            nvec[free] = vals
            if int(d @ nvec) != N: continue
            if (A @ nvec >= 0).all() and (nvec @ nvec) <= N * N:
                found.append(nvec.copy())
        if found:
            return N, found
    return None, []

for target in ['2', '3', '4', "4'", '6', '5']:
    N, found = minimal_configs(target, Nmax=4)
    print(f"\n=== pole of L({target}) : minimal ord zeta_L = {N}; {len(found)} configuration(s) ===")
    for f in found:
        print("   ", {k: int(v) for k, v in zip(names, f) if v != 0})
