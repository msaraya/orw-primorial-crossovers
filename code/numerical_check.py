# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Mohamed Hassan Ali Awaad Saraya

"""Reproduce the numerical consistency table in the paper.

Requires: mpmath
The computation is independent of the formal coefficient derivation.
"""
import mpmath as mp

mp.mp.dps = 80


def primes_below(q):
    ps = []
    for n in range(2, q):
        isprime = True
        for p in ps:
            if p*p > n:
                break
            if n % p == 0:
                isprime = False
                break
        if isprime:
            ps.append(n)
    return ps


def primorial_before(q):
    P = 1
    for p in primes_below(q):
        P *= p
    return P


def log_normalized_endpoint(D, L):
    """log J_D(L), where I_D(e^L)=e^L e^{-D/L} L^{-2} J_D(L)."""
    x = mp.mpf(D) / L
    upper = L - mp.log(2)

    def integrand(s):
        den = 1 - s/L
        return mp.e**(-s) * den**(-2) * mp.e**(-x*(1/den - 1))

    # The e^{-s} tail beyond 100 is negligible at 80-digit working precision
    # for the displayed digits. Increase cutoff to verify stability if desired.
    cutoff = min(upper, mp.mpf(100))
    nodes = [mp.mpf(0), 1, 3, 7, 15, 30, 60, cutoff]
    nodes = sorted(set(x for x in nodes if x <= cutoff))
    J = mp.quad(integrand, nodes)
    return mp.log(J)


def F(L, P, q):
    return (mp.log(mp.mpf(q-1)/(q-2))
            - mp.mpf((q-1)*P)/L
            + log_normalized_endpoint(q*P, L)
            - log_normalized_endpoint(P, L))


def root_and_expansions(q):
    P = primorial_before(q)
    s = mp.log(mp.mpf(q-1)/(q-2))
    L0 = mp.mpf((q-1)*P) / s
    # The paper proves L0 < L_int < L0+2. A secant start inside this interval
    # is therefore deterministic and does not rely on the asymptotic series.
    Lint = mp.findroot(lambda L: F(L, P, q), (L0+mp.mpf('0.8'), L0+mp.mpf('1.2')))
    Delta = Lint - L0

    S = mp.mpf((q+1)*P)
    Q = mp.mpf(P)**2 * (q*q + q + 1)
    a1 = 1 + 3/L0
    a2 = a1 + (12-S/2)/L0**2
    a4 = (a2
          + (63-4*S)/L0**3
          + (420-mp.mpf(63)/2*S+Q/3)/L0**4)
    return P, L0, Delta, a1, a2, a4


def fmt(x):
    return mp.nstr(x, 11)


if __name__ == '__main__':
    print(' q      numerical Δ       1+3/L0        through L0^-2    through L0^-4')
    for q in (5, 7, 11, 13, 17):
        P, L0, Delta, a1, a2, a4 = root_and_expansions(q)
        print(f'{q:2d}  {fmt(Delta):>15}  {fmt(a1):>15}  {fmt(a2):>15}  {fmt(a4):>15}')
