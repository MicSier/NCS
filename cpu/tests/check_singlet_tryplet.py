#!/usr/bin/env python3
"""Independent check of singlet / tryplet (Tc, mu) points produced by NCS.

Re-solves the gap + occupation equations with scipy (adaptive quad with the
Fermi points as breakpoints, Newton/hybrid root-finding) - no code shared with
NCS.h. Usage:
    python3 check_singlet_tryplet.py singlet singlet.txt [Vs=2 nl=1.1 tp=0]
    python3 check_singlet_tryplet.py tryplet tryplet.txt [Vt=2 nl=1.1 tp=0]
Each line of the .txt is: gamma0 Tc mu
"""
import sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, fsolve


def xi(x, mu, tp):
    return -2.0 * np.cos(x) + 2.0 * tp * np.cos(2.0 * x) - mu


def g(x, g0):
    return g0 * abs(np.sin(x))


def h(e, T):  # tanh(e/2T)/(2e), even, h(0) = 1/(4T)
    return 1.0 / (4.0 * T) if abs(e) < 1e-12 else np.tanh(e / (2.0 * T)) / (2.0 * e)


def fermi_points(mu, g0, tp):
    xs = np.linspace(0.0, np.pi, 4001)
    pts = []
    for s in (+1.0, -1.0):
        f = lambda x: xi(x, mu, tp) + s * g(x, g0)
        v = np.array([f(x) for x in xs])
        for i in np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]:
            pts.append(brentq(f, xs[i], xs[i + 1], xtol=1e-15))
    return sorted(p for p in pts if 0.0 < p < np.pi)


def integ(f, pts):
    val, _ = quad(f, 0.0, np.pi, points=pts or None, limit=2000, epsabs=1e-14, epsrel=1e-13)
    return val


def occupation(T, mu, g0, tp, nl, pts):
    f = lambda x: np.tanh((xi(x, mu, tp) + g(x, g0)) / (2 * T)) + np.tanh((xi(x, mu, tp) - g(x, g0)) / (2 * T))
    return 1.0 - nl - integ(f, pts) / (2.0 * np.pi)


def singlet_gap(T, mu, g0, tp, V, pts):
    f = lambda x: h(xi(x, mu, tp) + g(x, g0), T) + h(xi(x, mu, tp) - g(x, g0), T)
    return 1.0 - V / (2.0 * np.pi) * integ(f, pts)


def tryplet_gap(T, mu, g0, tp, V, pts):
    # d(k) || x-direction with |d|^2 = sin^2 k; gamma(k) = g0 sin k (projected on d)
    def f(x):
        e, vg, d2 = xi(x, mu, tp), g(x, g0), np.sin(x) ** 2
        hp, hm = h(e + vg, T), h(e - vg, T)
        u1 = hp + hm
        if vg < 1e-12:
            return u1 * d2
        gd = g0 * np.sin(x) ** 2
        pare = d2 - (gd / vg) ** 2
        u2 = (vg / e) * (hp - hm) if abs(e) > 1e-9 else 2 * vg * (
            (1 - np.tanh(vg / (2 * T)) ** 2) / (4 * T * vg) - np.tanh(vg / (2 * T)) / (2 * vg * vg))
        return u1 * d2 + u2 * pare
    return 1.0 - V / (2.0 * np.pi) * integ(f, pts)


def solve(kind, g0, T0, mu0, V, nl, tp):
    gap = singlet_gap if kind == "singlet" else tryplet_gap

    def F(z):
        T, mu = z
        pts = fermi_points(mu, g0, tp)
        return [gap(T, mu, g0, tp, V, pts), occupation(T, mu, g0, tp, nl, pts)]

    sol, info, ier, msg = fsolve(F, [T0, mu0], full_output=True, xtol=1e-13)
    return sol, ier == 1, np.max(np.abs(info["fvec"]))


def main():
    kind, path = sys.argv[1], sys.argv[2]
    kv = dict(a.split("=") for a in sys.argv[3:])
    V = float(kv.get("Vs", kv.get("Vt", 2.0)))
    nl, tp = float(kv.get("nl", 1.1)), float(kv.get("tp", 0.0))
    rows = np.loadtxt(path, ndmin=2)
    step = max(1, len(rows) // 8)
    sel = list(rows[::step]) + ([rows[-1]] if len(rows) % step != 1 else [])
    print(f"{'gamma0':>12} {'Tc (NCS)':>14} {'Tc (scipy)':>14} {'rel dTc':>9} {'mu (NCS)':>12} {'mu (scipy)':>12} {'|d mu|':>9}")
    worst = 0.0
    for g0, T, mu in sel:
        (Ts, mus), ok, res = solve(kind, g0, T, mu, V, nl, tp)
        rel = abs(Ts - T) / T
        worst = max(worst, rel)
        print(f"{g0:12.8f} {T:14.10f} {Ts:14.10f} {rel:9.1e} {mu:12.8f} {mus:12.8f} {abs(mus - mu):9.1e}"
              + ("" if ok else f"   (fsolve not converged, |F|={res:.1e})"))
    print(f"worst relative Tc deviation: {worst:.2e}")


if __name__ == "__main__":
    main()
