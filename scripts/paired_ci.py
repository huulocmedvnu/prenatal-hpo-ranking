"""Paired difference of two proportions with a 95% confidence interval.
b = fetuses in the top k only with the second method, c = only with the first, n = fetuses. Difference = (b - c) / n.
Tango (1998) score interval, checked against R PropCIs::scoreci.mp."""
import math
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import norm


def tango(b, c, n, level=0.95):
    """Score interval for p12 - p21 in a paired 2x2 table (Tango 1998), with the nuisance parameter p21
    estimated by constrained maximum likelihood."""
    z = norm.ppf(0.5 + level / 2)
    a = n - b - c  # concordant cells pooled, the likelihood depends on them only through their sum

    def p21_hat(delta):
        lo, hi = max(0.0, -delta), (1 - delta) / 2
        if hi - lo < 1e-12:
            return lo
        def nll(q):
            q12, q21 = q + delta, q
            qa = 1 - q12 - q21
            out = 0.0
            for k, p in ((b, q12), (c, q21), (a, qa)):
                if k:
                    if p <= 0:
                        return 1e300
                    out -= k * math.log(p)
            return out
        return minimize_scalar(nll, bounds=(lo, hi), method="bounded", options={"xatol": 1e-12}).x

    def stat(delta):
        q = p21_hat(delta)
        v = n * (2 * q + delta - delta * delta)
        return (b - c - n * delta) / math.sqrt(max(v, 1e-300))
    d = (b - c) / n
    eps = 1e-9
    lo = -1 + eps if stat(-1 + eps) < z else brentq(lambda x: stat(x) - z, -1 + eps, d)
    hi = 1 - eps if stat(1 - eps) > -z else brentq(lambda x: stat(x) + z, d, 1 - eps)
    return d, lo, hi


def paired(t1, t2, level=0.95):
    """t1, t2: boolean sequences (first and second method). Returns counts and both intervals in percentage points."""
    t1, t2 = list(map(bool, t1)), list(map(bool, t2))
    n = len(t1); x1, x2 = sum(t1), sum(t2)
    b = sum(y and not x for x, y in zip(t1, t2)); c = sum(x and not y for x, y in zip(t1, t2))
    d, tl, tu = tango(b, c, n, level)
    r = lambda v: round(100 * v, 1)
    return dict(n=n, first=x1, second=x2, gained=b, lost=c, diff_pp=r(d), ci=[r(tl), r(tu)])


if __name__ == "__main__":
    # Reference values from R PropCIs::scoreci.mp(6, 5, 229) = (-0.03723528, 0.02758603) for (c - b) / n
    d, lo, hi = tango(6, 5, 229)
    assert abs(lo + 0.02758603) < 1e-6 and abs(hi - 0.03723528) < 1e-6, (lo, hi)
    print("ok", d, lo, hi)
