"""
This file checks prime
powers q=p^s with p>=5, 17 <= q < 53267, and

    (1 + sqrt(8q - 7)) / 2 < k <= q / 2.

The corrected large-range inequality is

    E(k,q,p) = (q)_k - (q-1)L_{k,q,p} - (q^2-q)Q_{k,q,p} > 0,

where

    Q_{k,q,p} = (sqrt(q) + k + (q - sqrt(q))/p - 1)_k,
    L_{k,q,p} = k! * binom(q/p + k/p - 1, k/p) if p divides k,
                0 otherwise.
"""
from __future__ import annotations

import argparse
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from decimal import Decimal, getcontext
from os import cpu_count
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

Q_MIN = 17
Q_MAX_EXCLUSIVE = 53267
CHAR_MIN = 5

Pair = Tuple[int, int]


@dataclass(frozen=True)
class CaseResult:
    q: int
    p: int
    s: int
    status: str
    k_min: Optional[int]
    k_max: Optional[int]
    total_k: int
    hold_count: int
    fail_count: int
    fail_ranges: Tuple[Tuple[int, int], ...]


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    if n % 3 == 0:
        return n == 3
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def generate_prime_powers(lo: int, hi_exclusive: int, char_min: int) -> List[Tuple[int, int, int]]:
    cases: List[Tuple[int, int, int]] = []
    seen: Set[int] = set()
    for p in range(char_min, hi_exclusive):
        if not is_prime(p):
            continue
        q, s = p, 1
        while q < hi_exclusive:
            if q >= lo and q not in seen:
                seen.add(q)
                cases.append((q, p, s))
            q *= p
            s += 1
    cases.sort()
    return cases


def k_lower(q: int) -> int:
    """Smallest integer k with (1+sqrt(8q-7))/2 < k, exactly."""
    return (1 + math.isqrt(8 * q - 7)) // 2 + 1


def k_upper(q: int) -> int:
    return q // 2


def logaddexp(x: float, y: float) -> float:
    if x == float("-inf"):
        return y
    if y == float("-inf"):
        return x
    m = max(x, y)
    return m + math.log(math.exp(x - m) + math.exp(y - m))


def log_binom(n: int, r: int) -> float:
    if r < 0 or r > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(r + 1) - math.lgamma(n - r + 1)


def log_L_term(q: int, p: int, k: int) -> float:
    if k % p != 0:
        return float("-inf")
    r = k // p
    n = q // p + r - 1
    return math.lgamma(k + 1) + log_binom(n, r)


def initial_log_falling_q(q: int, k: int) -> float:
    return sum(math.log(q - i) for i in range(k))


def initial_log_Q(q: int, p: int, k: int) -> float:
    sqrt_q = math.sqrt(q)
    delta = (q - sqrt_q) / p
    a_k = sqrt_q + k + delta - 1.0
    return sum(math.log(a_k - i) for i in range(k))


def check_inequality_range(q: int, p: int, s: int, verbose: bool = False) -> Dict[int, bool]:
    """Return {k: E(k,q,p)>0} for the large range of a fixed q."""
    sqrt_q = math.sqrt(q)
    delta = (q - sqrt_q) / p
    k_min, k_max = k_lower(q), k_upper(q)
    if k_min > k_max:
        return {}

    A = initial_log_falling_q(q, k_min)  # log((q)_k)
    B = initial_log_Q(q, p, k_min)       # log(Q_{k,q,p})
    log_q_minus_1 = math.log(q - 1)
    log_q2_minus_q = math.log(q * (q - 1))

    out: Dict[int, bool] = {}
    for k in range(k_min, k_max + 1):
        log_L = log_L_term(q, p, k)
        C = float("-inf") if log_L == float("-inf") else log_q_minus_1 + log_L
        D = log_q2_minus_q + B
        rhs_log = logaddexp(C, D)
        out[k] = A > rhs_log
        if verbose and (not out[k] or k <= k_min + 5):
            print(f"q={q}={p}^{s}, k={k}, A={A:.12g}, C={C:.12g}, D={D:.12g}, holds={out[k]}")
        if k < k_max:
            A += math.log(q - k)
            B += math.log(sqrt_q + k + delta)
    return out


def collapse_to_ranges(values: Iterable[int]) -> Tuple[Tuple[int, int], ...]:
    vals = sorted(values)
    if not vals:
        return tuple()
    ranges: List[Tuple[int, int]] = []
    start = end = vals[0]
    for x in vals[1:]:
        if x == end + 1:
            end = x
        else:
            ranges.append((start, end))
            start = end = x
    ranges.append((start, end))
    return tuple(ranges)


def check_one_case(case: Tuple[int, int, int]) -> CaseResult:
    q, p, s = case
    k_min, k_max = k_lower(q), k_upper(q)
    if k_min > k_max:
        return CaseResult(q, p, s, "empty", None, None, 0, 0, 0, tuple())
    table = check_inequality_range(q, p, s)
    fail_ks = [k for k, ok in table.items() if not ok]
    total = len(table)
    fail_count = len(fail_ks)
    hold_count = total - fail_count
    status = "all_hold" if fail_count == 0 else ("all_fail" if hold_count == 0 else "some_fail")
    return CaseResult(q, p, s, status, k_min, k_max, total, hold_count, fail_count, collapse_to_ranges(fail_ks))


# Optional high-precision checks for selected triples.
def decimal_falling(x: Decimal, n: int) -> Decimal:
    out = Decimal(1)
    for i in range(n):
        out *= x - Decimal(i)
    return out


def decimal_direct_check(q: int, p: int, k: int, prec: int = 160) -> bool:
    getcontext().prec = prec
    qd, pd, kd = Decimal(q), Decimal(p), Decimal(k)
    sqrt_qd = qd.sqrt()
    lhs = decimal_falling(qd, k)
    Q_inner = sqrt_qd + kd + (qd - sqrt_qd) / pd - Decimal(1)
    rhs = (qd * qd - qd) * decimal_falling(Q_inner, k)
    if k % p == 0:
        r = k // p
        L = math.factorial(k) * math.comb(q // p + r - 1, r)
        rhs += Decimal(q - 1) * Decimal(L)
    return lhs > rhs


def verify_selected_points(results: Sequence[CaseResult]) -> int:
    mismatches = 0
    checked = 0
    for r in results:
        if r.status == "empty":
            continue
        assert r.k_min is not None and r.k_max is not None
        recurrence = check_inequality_range(r.q, r.p, r.s)
        ks = {r.k_min, r.k_max, (r.k_min + r.k_max) // 2}
        for a, b in r.fail_ranges:
            for x in (a - 1, a, b, b + 1):
                if r.k_min <= x <= r.k_max:
                    ks.add(x)
        for k in sorted(ks):
            expected = recurrence[k]
            actual = decimal_direct_check(r.q, r.p, k)
            checked += 1
            if expected != actual:
                print(f"MISMATCH: q={r.q}={r.p}^{r.s}, k={k}, recurrence={expected}, decimal={actual}")
                mismatches += 1
    print(f"Decimal verification checked {checked} selected k-values.")
    return mismatches


# Direct DP verification of finite exceptional cases.
def prime_field_tables(q: int) -> Tuple[List[List[int]], List[int]]:
    add = [[(i + j) % q for j in range(q)] for i in range(q)]
    square = [(i * i) % q for i in range(q)]
    return add, square


def gf25_tables() -> Tuple[List[List[int]], List[int]]:
    # F_25 = F_5[a]/(a^2-2), encoded as a0 + 5*a1.
    p, nonsquare, q = 5, 2, 25
    def enc(a: int, b: int) -> int:
        return (a % p) + p * (b % p)
    def dec(x: int) -> Tuple[int, int]:
        return x % p, x // p
    add = [[0] * q for _ in range(q)]
    square = [0] * q
    for x in range(q):
        a, b = dec(x)
        for y in range(q):
            c, d = dec(y)
            add[x][y] = enc(a + c, b + d)
        square[x] = enc(a * a + nonsquare * b * b, 2 * a * b)
    return add, square


def coverage_sizes(q: int, kmax: int, add: List[List[int]], square: List[int]) -> List[int]:
    states: List[Set[Pair]] = [set() for _ in range(kmax + 1)]
    states[0].add((0, 0))
    for idx, x in enumerate(range(q)):
        x2 = square[x]
        for r in range(min(idx, kmax - 1), -1, -1):
            for s1, s2 in tuple(states[r]):
                states[r + 1].add((add[s1][x], add[s2][x2]))
    return [len(layer) for layer in states]


def verify_exception_cases(failures: Set[Tuple[int, int]]) -> bool:
    """
    Verify exactly the (q,k) pairs that were found to fail the inequality.
    Uses direct DP to check that they cover all q^2 pairs (sum x, sum x^2).
    """
    if not failures:
        print("\nNo runtime failures to verify with DP.")
        return True

    by_q: Dict[int, List[int]] = {}
    for q, k in failures:
        by_q.setdefault(q, []).append(k)

    print("\n" + "=" * 80)
    print("DIRECT DP VERIFICATION OF RUNTIME FAILURES")
    print("=" * 80)
    ok_all = True

    for q, ks in by_q.items():
        if q == 25:
            add, square = gf25_tables()
        else:
            add, square = prime_field_tables(q)

        max_k = max(ks)
        sizes = coverage_sizes(q, max_k, add, square)

        for k in ks:
            ok = sizes[k] == q * q
            ok_all = ok_all and ok
            status = "OK" if ok else "FAIL"
            print(f"q={q:2d}, k={k:2d}: reached {sizes[k]:4d}/{q*q:4d} pairs -- {status}")
    return ok_all


def format_ranges(ranges: Tuple[Tuple[int, int], ...]) -> str:
    if not ranges:
        return "-"
    return ", ".join(f"{{{a}}}" if a == b else f"[{a}, {b}]" for a, b in ranges)


def print_report(results: Sequence[CaseResult], failures: Set[Tuple[int, int]]) -> None:
    all_hold = [r for r in results if r.status == "all_hold"]
    some_fail = [r for r in results if r.status == "some_fail"]
    all_fail = [r for r in results if r.status == "all_fail"]
    empty = [r for r in results if r.status == "empty"]

    print("\n" + "=" * 80)
    print("RESULTS FOR THE CORRECTED INEQUALITY (RUNTIME COMPUTATION)")
    print("=" * 80)
    print(f"Total prime-power q: {len(results)}")
    print(f"All hold:            {len(all_hold)}")
    print(f"Some fail:           {len(some_fail)}")
    print(f"All fail:            {len(all_fail)}")
    print(f"Empty k-range:       {len(empty)}")

    print("\n" + "=" * 80)
    print(f"RUNTIME FAILURE CASES  (total {len(failures)} points)")
    print("=" * 80)
    failures_by_case = sorted(some_fail + all_fail, key=lambda r: r.q)
    if not failures_by_case:
        print("None.")
    else:
        for r in failures_by_case:
            print(f"q={r.q:>5} = {r.p}^{r.s}, k in [{r.k_min}, {r.k_max}], "
                  f"holds {r.hold_count}/{r.total_k}, fails at {format_ranges(r.fail_ranges)}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Check the corrected large-k inequality.")
    parser.add_argument("--qmin", type=int, default=Q_MIN)
    parser.add_argument("--qmax-exclusive", type=int, default=Q_MAX_EXCLUSIVE)
    parser.add_argument("--char-min", type=int, default=CHAR_MIN)
    parser.add_argument("--workers", type=int, default=cpu_count() or 1)
    parser.add_argument("--decimal", action="store_true", help="run selected high-precision checks")
    parser.add_argument("--verify-exceptions", action="store_true",
                        help="directly verify the runtime failure cases by dynamic programming")
    parser.add_argument("--verbose", action="store_true", help="print detailed data for each checked k")
    args = parser.parse_args(argv)

    cases = generate_prime_powers(args.qmin, args.qmax_exclusive, args.char_min)
    print("=" * 80)
    print("CHECKING CORRECTED INEQUALITY")
    print("=" * 80)
    print(f"q range:            [{args.qmin}, {args.qmax_exclusive})")
    print(f"char(q) >=          {args.char_min}")
    print(f"Prime-power cases:  {len(cases)}")
    print(f"Workers:            {args.workers}")

    if args.verbose or args.workers <= 1:
        results = [check_one_case(case) for case in cases]
    else:
        results: List[CaseResult] = []
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(check_one_case, case) for case in cases]
            for i, future in enumerate(as_completed(futures), start=1):
                results.append(future.result())
                if i % 500 == 0 or i == len(futures):
                    print(f"  completed {i}/{len(futures)} cases")
    results.sort(key=lambda r: r.q)

    # Compute the runtime failure set
    failures = {(r.q, k) for r in results for a, b in r.fail_ranges for k in range(a, b + 1)}

    print_report(results, failures)

    if args.decimal:
        print("\n" + "=" * 80)
        print("SELECTED DECIMAL VERIFICATION")
        print("=" * 80)
        if verify_selected_points(results):
            return 1
        print("Selected Decimal verification passed.")

    if args.verify_exceptions:
        if not verify_exception_cases(failures):
            print("ERROR: Runtime failure cases did NOT pass DP verification!")
            return 1
        else:
            print("All runtime failure cases passed DP verification.")

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
