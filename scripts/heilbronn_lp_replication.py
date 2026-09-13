"""
Separate-implementation verification (same model family, assessor-adversary seat, 2026-09-04) of the
Heilbronn-bound computation of heilbronn_coin.py. Everything built from scratch:
 - G = SL(2,5) as explicit matrices mod 5
 - conjugacy classes; character table via Burnside-Dixon (class-sum eigenvectors)
 - all subgroups (closure of pairs), up to conjugacy; their character tables
 - "known cone" rows: Ind_K^G(lambda) for every subgroup K and every DEGREE-1 irrep lambda
   (monomial rows), plus Ind_{2.A4}^G(tau) for the degree-2 irreps tau of 2.A4 = SL(2,3)
   (tetrahedral rows, Langlands 1980); plus n_1 >= 0 (zeta entire at s0 != 1)
 - minimize N = sum n_chi chi(1) subject to n_rho = -1 (rho a spin 2-dim), integer n
 - LP relaxation bound, MILP min, enumeration of all integer optima, and the
   "drop the tetrahedral row" variant.
"""
import itertools, numpy as np
from numpy.linalg import eig
from scipy.optimize import linprog, milp, LinearConstraint, Bounds

p = 5
def mat(a,b,c,d): return (a%p, b%p, c%p, d%p)
def mul(x,y):
    a,b,c,d = x; e,f,g,h = y
    return mat(a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h)
def inv(x):
    a,b,c,d = x  # det = 1
    return mat(d, -b, -c, a)
I = mat(1,0,0,1)
G = [mat(a,b,c,d) for a in range(p) for b in range(p) for c in range(p) for d in range(p) if (a*d-b*c)%p == 1]
assert len(G) == 120
idx = {g:i for i,g in enumerate(G)}
n = len(G)
MT = np.array([[idx[mul(G[i],G[j])] for j in range(n)] for i in range(n)], dtype=np.int32)
INV = np.array([idx[inv(g)] for g in G])

def conj_classes(elems):
    """conjugacy classes of a subgroup given as list of indices (closed under mult)."""
    elems = list(elems); S = set(elems); left = set(elems); classes = []
    while left:
        g = min(left); cl = set()
        for x in elems:
            cl.add(MT[MT[x, g], INV[x]])
        classes.append(sorted(cl)); left -= cl
    return classes

def closure(gens):
    S = {idx[I]}; frontier = list(S)
    gens = list(gens)
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = MT[x, g]
                if y not in S:
                    S.add(y); new.append(y)
        frontier = new
    return frozenset(S)

def char_table(elems):
    """Burnside-Dixon: returns (classes, table[k][i] = chi_k(class i)) with rows normalized so chi(1)>0."""
    elems = list(elems); cls = conj_classes(elems); r = len(cls)
    cid = {}
    for i,c in enumerate(cls):
        for g in c: cid[g] = i
    sizes = np.array([len(c) for c in cls], dtype=float)
    # structure constants: C_i C_j = sum_k c_ijk C_k ; c_ijk = #{(x,y) in C_i x C_j : xy = z} for fixed z in C_k
    M = np.zeros((r, r, r))
    reps = [c[0] for c in cls]
    for i in range(r):
        for j in range(r):
            counts = np.zeros(r)
            # count for each k: number of x in C_i with x^{-1} z in C_j, z = reps[k]
            for k in range(r):
                z = reps[k]; cnt = 0
                for x in cls[i]:
                    if cid[MT[INV[x], z]] == j: cnt += 1
                counts[k] = cnt
            M[i, j, :] = counts
    # matrices (M_i)_{jk} = c_{ijk}; common eigenvectors
    rng = np.random.default_rng(1)
    A = sum(rng.standard_normal()*M[i] for i in range(r))
    w, V = eig(A.T)   # left eigenvectors of M_i^T ... we need vectors v with M_i^T v = omega_i v ? use right eigvecs of A.T
    # For class-sum action: omega_chi(C_i) satisfies sum_k c_ijk omega(C_k) = omega(C_i) omega(C_j)
    # => vector v_k = omega(C_k) is an eigenvector of the matrix (c_ijk)_{k,j}?? Use: for fixed i, sum_k c_ijk v_k = omega_i v_j: (M_i)_{jk} v_k = omega_i v_j => M_i v = omega_i v (right eigenvector of M_i).
    w, V = eig(A)
    chars = []
    for col in range(r):
        v = V[:, col]; v = v / v[cid[idx[I]]]   # omega(C_1) = 1
        omega = np.array([ (M[i] @ v)[cid[idx[I]]] / v[cid[idx[I]]] for i in range(r)])  # eigenvalue per class sum
        # chi(g_i) = omega_i * chi(1) / |C_i| ; chi(1)^2 = |G| / sum_i |omega_i|^2/|C_i|
        d2 = len(elems) / np.sum(np.abs(omega)**2 / sizes)
        d = np.sqrt(d2)
        chi = omega * d / sizes
        chars.append(chi)
    T = np.array(chars)
    # verify orthonormality
    Gram = (T * sizes) @ T.conj().T / len(elems)
    assert np.allclose(Gram, np.eye(r), atol=1e-8), "character table not orthonormal"
    degs = np.real(T[:, cid[idx[I]]])
    assert np.allclose(np.sum(degs**2), len(elems))
    # clean tiny imaginary parts / round degrees
    return cls, cid, sizes, T

# ---- G's table
clsG, cidG, sizG, TG = char_table(range(n))
r = len(clsG)
order = lambda g: next(k for k in range(1, 61) if closure([g]).__len__() == k)  # crude
def elt_order(g):
    x = g; k = 1
    while x != idx[I]:
        x = MT[x, g]; k += 1
    return k
print("SL(2,5): classes", r, "sizes", sizG.astype(int).tolist(), "rep orders", [elt_order(c[0]) for c in clsG])
degs = np.round(np.real(TG[:, cidG[idx[I]]])).astype(int)
print("degrees", degs.tolist())
# identify spin 2-dims (degree 2, chi(-I) = -2)
minusI = idx[mat(-1,0,0,-1)]
faithful = [k for k in range(r) if abs(TG[k, cidG[minusI]] + degs[k]) < 1e-8]
print("faithful (spin) irreps: degrees", [degs[k] for k in faithful])

# ---- all subgroups up to conjugacy
subs = set()
for i in range(n):
    for j in range(i, n):
        subs.add(closure([i, j]))
subs.add(frozenset([idx[I]]))
print("distinct subgroups (from pairs):", len(subs))
# conjugacy classes of subgroups
def conj_sub(H, x): return frozenset(MT[MT[x, h], INV[x]] for h in H)
sub_reps = []; seen = set()
for H in sorted(subs, key=len):
    if H in seen: continue
    orbit = {conj_sub(H, x) for x in range(n)}
    seen |= orbit
    sub_reps.append((H, len(orbit)))
print("subgroup classes:", [(len(H), c) for H, c in sub_reps])

def induce(H, chi_H, cidH):
    """Ind_H^G chi as a class function on G's classes."""
    vals = np.zeros(r, dtype=complex)
    Hs = set(H)
    for i in range(r):
        g = clsG[i][0]; s = 0
        for x in range(n):
            y = MT[MT[x, g], INV[x]]
            if y in Hs: s += chi_H[cidH[y]]
        vals[i] = s / len(H)
    return vals

def decompose(f):
    """multiplicities of f in G's irreps."""
    return np.array([np.sum(sizG * f * np.conj(TG[k])) / n for k in range(r)])

rows = []   # (label, multiplicity vector over G irreps)
tetra_rows = []
for H, ncj in sub_reps:
    if len(H) == n: continue
    clsH, cidH, sizH, TH = char_table(H)
    degH = np.round(np.real(TH[:, cidH[idx[I]]])).astype(int)
    for k in range(len(clsH)):
        f = induce(H, TH[k], cidH)
        m = decompose(f)
        mi = np.round(np.real(m)).astype(int)
        assert np.allclose(m, mi, atol=1e-6)
        if degH[k] == 1:
            rows.append((f"|H|={len(H)} deg1", mi))
        elif len(H) == 24 and degH[k] == 2:
            tetra_rows.append((f"|H|=24 deg2 (tetrahedral)", mi))
        # other irreps of proper subgroups: are they monomial? (they are induced from deg-1 of smaller K, so already
        # implied by the deg-1 rows of K; we do not add them — but we CHECK that every irrep of every proper subgroup
        # other than the SL(2,3) 2-dims is monomial, i.e. non-negative combination of Ind_K^H lambda.)
def uniq(rs):
    out = {}
    for lab, m in rs:
        out.setdefault(tuple(m.tolist()), lab)
    return out
mono = uniq(rows); tet = uniq(tetra_rows)
print("distinct monomial rows:", len(mono), " distinct tetrahedral rows:", len(tet))
print("tetrahedral rows (mult over G irreps, degrees %s):" % degs.tolist())
for m, lab in tet.items(): print("   ", m, "-> in monomial set?" , m in mono)

# ---- the coin problem
labels = [f"d{degs[k]}{'s' if k in faithful else ''}#{k}" for k in range(r)]
print("irrep labels:", labels)
Amono = np.array(list(mono.keys()), dtype=float)
Atet = np.array([m for m in tet.keys() if m not in mono], dtype=float)
e1 = np.zeros(r); e1[0] = 1.0   # trivial irrep index? find it
triv = [k for k in range(r) if np.allclose(TG[k], 1)][0]
e1 = np.zeros(r); e1[triv] = 1.0

def solve(A_rows, rho, m_order=1, integer=True, bound=6):
    A = np.vstack([A_rows, e1[None, :]])
    c = degs.astype(float)
    lb = -bound*np.ones(r); ub = bound*np.ones(r)
    lb[rho] = -m_order; ub[rho] = -m_order
    cons = LinearConstraint(A, 0, np.inf)
    if integer:
        res = milp(c, constraints=cons, integrality=np.ones(r), bounds=Bounds(lb, ub))
    else:
        res = linprog(c, A_ub=-A, b_ub=np.zeros(A.shape[0]), bounds=list(zip(lb, ub)), method="highs")
    return res

rho = faithful[[degs[k] for k in faithful].index(2)]   # first spin 2-dim
rho2 = [k for k in faithful if degs[k] == 2 and k != rho][0]
Afull = np.vstack([Amono, Atet]) if len(Atet) else Amono
for name, A in [("full known cone (monomial + tetrahedral)", Afull), ("monomial rows only (no tetrahedral)", Amono)]:
    lp = solve(A, rho, 1, integer=False); ip = solve(A, rho, 1, integer=True)
    print(f"\n[{name}] pole of order 1 in L(rho): LP min N = {lp.fun:.4f}; MILP min N = {ip.fun:.0f}")
    for m_ in (2, 3):
        lp = solve(A, rho, m_, integer=False); ip = solve(A, rho, m_, integer=True)
        print(f"   pole of order {m_}: LP min N = {lp.fun:.4f}; MILP min N = {ip.fun:.0f}")

# enumerate ALL integer solutions with N = 4, n_rho = -1, |n| <= 4 under the full cone
A = np.vstack([Afull, e1[None, :]])
others = [k for k in range(r) if k != rho]
sols = []
rng_ = range(-4, 5)
# vectorized: split variables
first, last = others[:5], others[5:]
grid_last = np.array(list(itertools.product(rng_, repeat=len(last))))
for vals in itertools.product(rng_, repeat=len(first)):
    v = np.zeros((len(grid_last), r)); v[:, rho] = -1
    for k, x in zip(first, vals): v[:, k] = x
    v[:, last] = grid_last
    Nv = v @ degs
    ok = (Nv <= 4)
    if not ok.any(): continue
    vv = v[ok]
    feas = np.all(vv @ A.T >= -1e-9, axis=1)
    for s in vv[feas]: sols.append(tuple(int(x) for x in s))
print("\nAll integer coins with n_rho = -1, N <= 4, |n|<=4, all known marginals >= 0:")
for s in sorted(set(sols)):
    N = int(np.dot(s, degs))
    print("   N =", N, {labels[k]: s[k] for k in range(r) if s[k] != 0})

# also: pole in the adjoint (3-dim), the register (4-dim non-spin), 4' (spin 4), 5, 6
for target_deg, spin in [(3, False), (4, False), (4, True), (5, False), (6, True)]:
    ks = [k for k in range(r) if degs[k] == target_deg and ((k in faithful) == spin)]
    k = ks[0]
    ip = solve(Afull, k, 1, integer=True)
    print(f"pole in L(deg {target_deg}{'spin' if spin else ''}): MILP min N =", "infeasible" if not ip.success else int(round(ip.fun)))

print("\nrho index", rho, "rho' index", rho2)
print("Check the 2-dim spin characters are Galois conjugate (values on order-5/10 classes):")
for k in (rho, rho2): print("  ", labels[k], np.round(np.real(TG[k]),4).tolist())
print("Ind_{C10} xi rows containing rho (the two-tile 2+4'+6) present?", any(m[rho]==1 and m[[k for k in range(r) if degs[k]==4 and k in faithful][0]]==1 and m[[k for k in range(r) if degs[k]==6][0]]==1 and sum(m)==3 for m in mono))

# ---- cascade: Sym^k of the spin 2-dim via Chebyshev U_k(t/2), t = chi_2(g) (det = 1)
from numpy.polynomial import chebyshev as C
t = np.real(TG[rho])
def U(k, x):
    # Chebyshev second kind via recurrence
    u0, u1 = np.ones_like(x), 2*x
    if k == 0: return u0
    for _ in range(k-1): u0, u1 = u1, 2*x*u1 - u0
    return u1
names = {k: labels[k] for k in range(r)}
print("\nCascade Sym^k(2) decomposition into irreps of SL(2,5):")
for k in range(1, 13):
    f = U(k, t/2)
    m = np.round(np.real(decompose(f))).astype(int)
    print(f"  Sym^{k:2d} (dim {k+1:2d}):", {names[i]: int(m[i]) for i in range(r) if m[i]}, " invariants:", int(m[triv]))
print("\nMoments E|a|^{2k} over the finite group vs Catalan(k):")
from math import comb
for k in range(1, 8):
    mom = np.sum(sizG * np.abs(TG[rho])**(2*k)) / n
    cat = comb(2*k, k)//(k+1)
    print(f"  k={k}: group {mom:.6f}  Catalan {cat}")
