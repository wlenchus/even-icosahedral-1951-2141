"""
Reduce PARI/GP's side-pairing matrices for Gamma_0(N) (output of mspolygon, see
generators/mspolygon_N.gp) to a list of distinct generators up to sign and inverse.

Each side pairing gamma = [[a,b],[c,d]] appears in the raw list together with its inverse
(the paired side), and PARI may return either sign. The four matrices
    gamma, -gamma, gamma^-1 = [[d,-b],[-c,a]], -gamma^-1
generate the same subgroup, so we keep one representative per such class: the
lexicographically smallest of the four tuples. The translation T = [[1,1],[0,1]] is kept
(it is the 'trivial' generator, tested as 1-periodicity).

Usage:  python3 dedup_generators.py N
Reads   ../generators/side_pairings_raw_N.txt   (one matrix per line: a b c d)
Writes  ../generators/generators_unique_N.csv   (a,b,c,d; header line)
        ../generators/gens_unique_N.npy         (int64 array, shape (k,4)) -- the input of the test scripts
Expected: 652 -> 327 for N = 1951; 716 -> 359 for N = 2141.
"""
import sys, os, csv
import numpy as np

N = int(sys.argv[1])
here = os.path.dirname(os.path.abspath(__file__))
raw = os.path.join(here, '..', 'generators', f'side_pairings_raw_{N}.txt')

def canon(m):
    a, b, c, d = m
    assert a * d - b * c == 1, m
    reps = [(a, b, c, d), (-a, -b, -c, -d), (d, -b, -c, a), (-d, b, c, -a)]
    return min(reps)

seen = {}
order = []
with open(raw) as f:
    for line in f:
        parts = line.split()
        if len(parts) != 4:
            continue
        m = tuple(int(x) for x in parts)
        key = canon(m)
        if key not in seen:
            seen[key] = m          # keep the first-encountered representative, as PARI emitted it
            order.append(key)

gens = [seen[k] for k in order]
# sanity: every generator lies in Gamma_0(N)
for a, b, c, d in gens:
    assert c % N == 0 and a * d - b * c == 1

out_csv = os.path.join(here, '..', 'generators', f'generators_unique_{N}.csv')
out_npy = os.path.join(here, '..', 'generators', f'gens_unique_{N}.npy')
with open(out_csv, 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['a', 'b', 'c', 'd'])
    for g in gens: w.writerow(g)
np.save(out_npy, np.array(gens, dtype=np.int64))
n_raw = sum(1 for line in open(raw) if len(line.split()) == 4)
print(f"N = {N}: {n_raw} raw side pairings -> {len(gens)} distinct generators up to sign and inverse")
