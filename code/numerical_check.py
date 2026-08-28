#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Mohamed Hassan Ali Awaad Saraya

"""Reproduce the numerical consistency table in the manuscript.

This script evaluates the normalized ORW endpoint integral with high-precision
mpmath quadrature and solves the crossover in the proved interval (L0,L0+2).
It is ancillary only; no proof depends on the computation.
"""
import mpmath as mp

mp.mp.dps = 80


def primes_below(q):
    out=[]
    for n in range(2,q):
        ok=True
        for p in out:
            if p*p>n:
                break
            if n%p==0:
                ok=False
                break
        if ok:
            out.append(n)
    return out


def primorial_before(q):
    P=1
    for p in primes_below(q):
        P*=p
    return P


def log_J(D,L):
    """log J_D(1/L), where I_D(e^L)=e^(L-D/L)L^-2 J_D(1/L)."""
    D=mp.mpf(D); L=mp.mpf(L)
    x=D/L
    upper=L-mp.log(2)
    def f(s):
        den=1-s/L
        return mp.e**(-s) * den**(-2) * mp.e**(-x*s/L/den)
    # The omitted tail beyond 120 is far below the displayed precision.
    cutoff=min(upper, mp.mpf(120))
    pts=[mp.mpf(0),1,3,7,15,30,60,90,cutoff]
    pts=sorted(set(p for p in pts if p<=cutoff))
    J=mp.quad(f, pts)
    return mp.log(J)


def F(L,P,q):
    A=(q-1)*P
    s=mp.log(mp.mpf(q-1)/(q-2))
    return s-mp.mpf(A)/L + log_J(q*P,L)-log_J(P,L)


def row(q):
    P=primorial_before(q)
    s=mp.log(mp.mpf(q-1)/(q-2))
    L0=mp.mpf((q-1)*P)/s
    root=mp.findroot(lambda L:F(L,P,q), (L0+mp.mpf('0.5'), L0+mp.mpf('1.5')))
    delta=root-L0
    S=mp.mpf((q+1)*P)
    Q=mp.mpf(P)**2*(q*q+q+1)
    t1=1+3/L0
    t2=t1+(12-S/2)/L0**2
    t4=t2+(63-4*S)/L0**3+(420-mp.mpf(63)/2*S+Q/3)/L0**4
    return delta,t1,t2,t4


def fmt(x):
    return mp.nstr(x, 11)


if __name__=='__main__':
    print('q   numerical Delta   1+3/L0        through L0^-2   through L0^-4')
    for q in (5,7,11,13,17):
        vals=row(q)
        print(f'{q:<2}  ' + '  '.join(f'{fmt(v):>13}' for v in vals))
