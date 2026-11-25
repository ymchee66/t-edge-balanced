"""Forced counts for a 3-edge-balanced graph in G(n,k) (Section 3)."""
from fractions import Fraction as Fr
def targets(n, k):
    a = Fr(2*k*(k-1), n+1)
    c = Fr(4*k*(k-1)*(k-2), 3*(n+1)*(n*n-n-4))
    p = 3*(n-3)*c; s = (n-3)*c
    return dict(a=a, c=c, p=p, s=s, R=p+3*c)
