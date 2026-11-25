#!/usr/bin/env python3
"""Randomized local search for a graphical degree sequence d_1..d_n with
sum d_i = 2k, sum C(d_i,2) = a*, sum C(d_i,3) = s*.
Usage: degseq.py n k seed"""
import sys, random
from math import comb
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from targets import targets

def erdos_gallai(d):
    d = sorted(d, reverse=True); n = len(d)
    if sum(d) % 2: return False
    pref = 0
    for r in range(1, n + 1):
        pref += d[r-1]
        if pref > r*(r-1) + sum(min(x, r) for x in d[r:]): return False
    return True

n, k, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
rng = random.Random(seed)
T = targets(n, k)
A, S = int(T['a']), int(T['s'])
mu = 2*k/n
var = (2*A + 2*k)/n - mu*mu
d = [min(n-1, max(0, round(rng.gauss(mu, var**0.5)))) for _ in range(n)]
diff = 2*k - sum(d)
while diff:
    i = rng.randrange(n)
    if diff > 0 and d[i] < n-1: d[i] += 1; diff -= 1
    elif diff < 0 and d[i] > 0: d[i] -= 1; diff += 1
a = sum(comb(x, 2) for x in d); s = sum(comb(x, 3) for x in d)
def err(a, s): return 3*mu*abs(a - A) + abs(s - S)
e = err(a, s); it = 0

def moves(cnt):
    """Map (change in a, change in s) -> (x, y) for moving one unit of degree
    from a vertex of degree y to a vertex of degree x."""
    vals = [v for v in cnt if cnt[v] > 0]
    M = {}
    for x in vals:
        if x >= n-1: continue
        for y in vals:
            if y <= 0 or (x == y and cnt[x] < 2): continue
            M.setdefault((x - y + 1, comb(x, 2) - comb(y-1, 2)), (x, y))
    return M

def apply(cnt, x, y):
    for v, dv in ((x, -1), (x+1, 1), (y, -1), (y-1, 1)):
        cnt[v] = cnt.get(v, 0) + dv

def finish(d, a, s):
    """Exact finishing: try one move, then all pairs of moves, reaching (A, S)."""
    from collections import Counter
    cnt = Counter(d); need = (A - a, S - s)
    M = moves(cnt)
    if need in M: return [M[need]]
    for (da, ds), (x, y) in M.items():
        c2 = Counter(cnt); apply(c2, x, y)
        rest = (need[0] - da, need[1] - ds)
        M2 = None
        if rest in M:
            x2, y2 = M[rest]
            if c2.get(x2, 0) > 0 and c2.get(y2, 0) > 0 and x2 < n-1 and y2 > 0 and (x2 != y2 or c2[x2] >= 2):
                return [(x, y), (x2, y2)]
    return None

def realize(d, plan):
    for x, y in plan:
        i = next(t for t in range(n) if d[t] == x)
        j = next(t for t in range(n) if d[t] == y and t != i)
        d[i] += 1; d[j] -= 1

while e > 0:
    it += 1
    if it > 50_000_000: sys.exit("no degree sequence found")
    if it % 20000 == 0:
        plan = finish(d, a, s)
        if plan:
            realize(d, plan)
            a = sum(comb(x, 2) for x in d); s = sum(comb(x, 3) for x in d); e = err(a, s)
            break
    i, j = rng.randrange(n), rng.randrange(n)
    if i == j or d[i] >= n-1 or d[j] <= 0: continue
    na = a + d[i] - (d[j] - 1)                      # C(x+1,2)-C(x,2)=x, C(y-1,2)-C(y,2)=-(y-1)
    ns = s + comb(d[i], 2) - comb(d[j]-1, 2)        # C(x+1,3)-C(x,3)=C(x,2), C(y-1,3)-C(y,3)=-C(y-1,2)
    ne = err(na, ns)
    if ne <= e or rng.random() < 0.01:
        d[i] += 1; d[j] -= 1; a, s, e = na, ns, ne
assert sum(d) == 2*k and a == A and s == S and erdos_gallai(d)
print(' '.join(map(str, sorted(d, reverse=True))))
