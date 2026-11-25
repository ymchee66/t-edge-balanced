#!/usr/bin/env python3
"""Admissible parameter sets (n,k): n >= 6, 4 <= k <= C(n,2)/2, (C1) and (C2),
and all forced counts of graphs with two and three edges integral.
Usage: python3 admissible.py LO HI   (requires sympy)"""
import sys
from math import comb, perm, gcd
from itertools import product
from sympy import factorint
T={2:[(2,3),(8,4)],3:[(6,3),(2,4),(6,4),(4,5),(48,6)]}
def ok(n,k,t):
    N=comb(n,2); ck=comb(k,t); cN=comb(N,t)
    return all((ck*perm(n,v))%(cN*a)==0 for a,v in T[t] if v<=n)
lo,hi=int(sys.argv[1]),int(sys.argv[2])
for n in range(max(lo,6),hi+1):
    N=comb(n,2); C=(n+1)*(n*n-n-4)
    while C%2==0: C//=2
    while C%3==0: C//=3
    mods=[p**e for p,e in factorint(C).items()]
    ks=set()
    for ch in product(range(3),repeat=len(mods)):
        r,M=0,1
        for m,c in zip(mods,ch):
            r+=M*(((c-r)*pow(M,-1,m))%m); M*=m
        for k in range(r,N//2+1,C):
            if k>=4 and (2*k*(k-1))%(n+1)==0 and ok(n,k,2) and ok(n,k,3): ks.add(k)
    for k in sorted(ks): print(n,k,flush=True)
