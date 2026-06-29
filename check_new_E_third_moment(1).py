"""
Check prime powers q=p^s with p>=5, 17 <= q < 53267, and

    (1 + sqrt(8q - 7)) / 2 < k <= q / 2.

This version checks the NEW third-moment lower-bound inequality

    E(k,p,q) > 0,

where

    E(k,p,q) = (q)_k/q^3
               - ((q^3-q)/q^3) * L_{k,q,p}
               - ((q^3-q)/q^3) * Q_{k,q,p}.

Equivalently, to avoid underflow/overflow, it checks

    (q)_k > (q^3-q) * (L_{k,q,p} + Q_{k,q,p}).

Here

    L_{k,q,p} = k! * binom(q/p + k/p - 1, k/p), if p divides k,
                0, otherwise,

and

    Q_{k,q,p} = (2*sqrt(q) + k + (q - 2*sqrt(q))/p - 1)_k.

For EXPECTED_EXCEPTIONAL_CASES, the direct finite check is the THIRD moment
subset-sum question: coverage of all triples

    (sum x, sum x^2, sum x^3) in F_q^3,

not the old second-moment coverage of pairs (sum x, sum x^2).

The third-moment DP uses numpy boolean arrays and finite-field arithmetic over
GF(p^s). It is intended for the small exceptional cases only.
"""
from __future__ import annotations

import argparse
import itertools
import math
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from decimal import Decimal, getcontext
from os import cpu_count
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

Q_MIN = 17
Q_MAX_EXCLUSIVE = 53267
CHAR_MIN = 5

Pair = Tuple[int, int]
RangeMap = Dict[int, Tuple[Tuple[int, int], ...]]


# These are the analytic failures for the new E(k,p,q) over the default range,
# grouped as q -> k-ranges. They are useful as a regression check.
EXPECTED_ANALYTIC_FAILURE_RANGES: RangeMap = {
    17: ((7, 8),),
    19: ((7, 9),),
    23: ((8, 11),),
    25: ((8, 12),),
    29: ((9, 14),),
    31: ((9, 15),),
    37: ((10, 18),),
    41: ((10, 20),),
    43: ((10, 21),),
    47: ((11, 23),),
    49: ((11, 24),),
    53: ((11, 26),),
    59: ((12, 14),),
    61: ((12, 13),),
    125: ((17, 17),),
}

# Direct third-moment DP shows that (q,k)=(17,7) does NOT cover all of F_17^3:
# it reaches 4896/4913 triples. Therefore it is intentionally not included in
# EXPECTED_EXCEPTIONAL_CASES below. If you want to verify every analytic failure
# by DP, run with --verify-runtime-failures; the script will report this failure.
KNOWN_THIRD_MOMENT_NONCOVERAGE: Set[Pair] = {(17, 7)}


def expand_ranges(ranges_by_q: RangeMap) -> Set[Pair]:
    return {
        (q, k)
        for q, ranges in ranges_by_q.items()
        for a, b in ranges
        for k in range(a, b + 1)
    }


EXPECTED_ANALYTIC_FAILURES: Set[Pair] = expand_ranges(EXPECTED_ANALYTIC_FAILURE_RANGES)

# These are the analytic failures that are expected to be rescued by the direct
# third-moment subset-sum DP check.
EXPECTED_EXCEPTIONAL_CASES: Set[Pair] = EXPECTED_ANALYTIC_FAILURES - KNOWN_THIRD_MOMENT_NONCOVERAGE


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


@dataclass(frozen=True)
class FieldData:
    q: int
    p: int
    s: int
    square: Tuple[int, ...]
    cube: Tuple[int, ...]
    modulus: Optional[Tuple[int, ...]]


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


def factor_prime_power(q: int) -> Tuple[int, int]:
    for p in range(2, q + 1):
        if not is_prime(p) or q % p != 0:
            continue
        t, s = p, 1
        while t < q:
            t *= p
            s += 1
        if t == q:
            return p, s
    raise ValueError(f"q={q} is not a prime power")


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


def q_term_base(q: int, p: int) -> float:
    """Base B such that Q_{k,q,p} = B(B+1)...(B+k-1)."""
    sqrt_q = math.sqrt(q)
    return 2.0 * sqrt_q + (q - 2.0 * sqrt_q) / p


def initial_log_Q(q: int, p: int, k: int) -> float:
    base = q_term_base(q, p)
    return sum(math.log(base + i) for i in range(k))


def check_inequality_range(q: int, p: int, s: int, verbose: bool = False) -> Dict[int, bool]:
    """Return {k: E(k,p,q)>0} for the large range of a fixed q."""
    k_min, k_max = k_lower(q), k_upper(q)
    if k_min > k_max:
        return {}

    base = q_term_base(q, p)
    A = initial_log_falling_q(q, k_min)  # log((q)_k)
    B = initial_log_Q(q, p, k_min)       # log(Q_{k,q,p})

    # This follows the E definition in the prompt: both L and Q are multiplied
    # by q^3-q after clearing the common q^3 denominator.
    log_coeff = math.log(q**3 - q)

    out: Dict[int, bool] = {}
    for k in range(k_min, k_max + 1):
        log_L = log_L_term(q, p, k)
        log_L_part = float("-inf") if log_L == float("-inf") else log_coeff + log_L
        log_Q_part = log_coeff + B
        rhs_log = logaddexp(log_L_part, log_Q_part)
        out[k] = A > rhs_log

        if verbose and (not out[k] or k <= k_min + 5):
            print(
                f"q={q}={p}^{s}, k={k}, "
                f"log_lhs={A:.12g}, log_L_part={log_L_part:.12g}, "
                f"log_Q_part={log_Q_part:.12g}, holds={out[k]}"
            )

        if k < k_max:
            A += math.log(q - k)
            B += math.log(base + k)
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


def decimal_direct_check(q: int, p: int, k: int, prec: int = 240) -> bool:
    getcontext().prec = prec
    qd, pd, kd = Decimal(q), Decimal(p), Decimal(k)
    sqrt_qd = qd.sqrt()

    lhs = decimal_falling(qd, k)

    q_inner = Decimal(2) * sqrt_qd + kd + (qd - Decimal(2) * sqrt_qd) / pd - Decimal(1)
    rhs = Decimal(q**3 - q) * decimal_falling(q_inner, k)

    if k % p == 0:
        r = k // p
        L = math.factorial(k) * math.comb(q // p + r - 1, r)
        rhs += Decimal(q**3 - q) * Decimal(L)

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


# Finite-field arithmetic for direct third-moment DP.
def decode_base_p(x: int, p: int, s: int) -> Tuple[int, ...]:
    coeffs: List[int] = []
    for _ in range(s):
        coeffs.append(x % p)
        x //= p
    return tuple(coeffs)


def encode_base_p(coeffs: Sequence[int], p: int) -> int:
    out = 0
    place = 1
    for c in coeffs:
        out += (c % p) * place
        place *= p
    return out


def poly_mod(poly: Sequence[int], divisor: Sequence[int], p: int) -> List[int]:
    """Return poly mod divisor over F_p. Coefficients are low-to-high."""
    a = [x % p for x in poly]
    d = len(divisor) - 1
    inv_lead = pow(divisor[-1], -1, p)
    for j in range(len(a) - 1, d - 1, -1):
        c = a[j] % p
        if c == 0:
            continue
        factor = (c * inv_lead) % p
        for i in range(d + 1):
            a[j - d + i] = (a[j - d + i] - factor * divisor[i]) % p
    out = a[:d]
    while out and out[-1] == 0:
        out.pop()
    return out


def is_irreducible_monic(poly: Sequence[int], p: int) -> bool:
    """Brute-force irreducibility test, adequate for small exceptional fields."""
    degree = len(poly) - 1
    if degree <= 0 or poly[-1] != 1:
        return False
    for d in range(1, degree // 2 + 1):
        for coeffs in itertools.product(range(p), repeat=d):
            divisor = list(coeffs) + [1]
            if not poly_mod(poly, divisor, p):
                return False
    return True


def find_irreducible_monic(p: int, s: int) -> Tuple[int, ...]:
    """Find a monic irreducible polynomial of degree s over F_p."""
    if s == 1:
        return (0, 1)
    for coeffs in itertools.product(range(p), repeat=s):
        # Constant term 0 would be divisible by x.
        if coeffs[0] == 0:
            continue
        candidate = tuple(coeffs) + (1,)
        if is_irreducible_monic(candidate, p):
            return candidate
    raise ValueError(f"Could not find an irreducible polynomial for GF({p}^{s})")


def build_field_data(q: int, p: Optional[int] = None, s: Optional[int] = None) -> FieldData:
    if p is None or s is None:
        p, s = factor_prime_power(q)
    if p**s != q:
        raise ValueError(f"Invalid field data: q={q}, p={p}, s={s}")

    if s == 1:
        square = tuple((x * x) % p for x in range(p))
        cube = tuple((x * x * x) % p for x in range(p))
        return FieldData(q=q, p=p, s=s, square=square, cube=cube, modulus=None)

    modulus = find_irreducible_monic(p, s)
    elements = [decode_base_p(x, p, s) for x in range(q)]

    def mul_idx(a_idx: int, b_idx: int) -> int:
        a = elements[a_idx]
        b = elements[b_idx]
        prod = [0] * (2 * s - 1)
        for i, ai in enumerate(a):
            if ai == 0:
                continue
            for j, bj in enumerate(b):
                prod[i + j] = (prod[i + j] + ai * bj) % p

        # Reduce by modulus x^s + modulus[s-1]x^{s-1}+...+modulus[0].
        for d in range(2 * s - 2, s - 1, -1):
            c = prod[d] % p
            if c == 0:
                continue
            for j in range(s):
                prod[d - s + j] = (prod[d - s + j] - c * modulus[j]) % p
        return encode_base_p(prod[:s], p)

    square = tuple(mul_idx(x, x) for x in range(q))
    cube = tuple(mul_idx(square[x], x) for x in range(q))
    return FieldData(q=q, p=p, s=s, square=square, cube=cube, modulus=modulus)


def coverage_sizes_moment3_numpy(
    field: FieldData,
    kmax: int,
    target_ks: Optional[Set[int]] = None,
    quiet: bool = False,
) -> List[int]:
    """
    Exact DP for third-moment subset sums.

    It computes the number of reachable triples in F_q^3 for each subset size
    0..kmax, where every selected field element is used at most once.
    """
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("numpy is required for third-moment DP verification") from exc

    q, p, s = field.q, field.p, field.s
    shape = (p,) * (3 * s)
    axes = tuple(range(3 * s))
    zero = (0,) * (3 * s)
    dec = [decode_base_p(x, p, s) for x in range(q)]

    states = [np.zeros(shape, dtype=np.bool_) for _ in range(kmax + 1)]
    states[0][zero] = True

    target_size = q**3
    watched = set(target_ks or [])
    for idx, x in enumerate(range(q)):
        shift = dec[x] + dec[field.square[x]] + dec[field.cube[x]]
        for r in range(min(idx, kmax - 1), -1, -1):
            states[r + 1] |= np.roll(states[r], shift=shift, axis=axes)

        if watched and all(int(states[k].sum()) == target_size for k in watched):
            if not quiet:
                print(f"  q={q}: all requested k-values became full after processing {idx + 1}/{q} elements")
            break

    return [int(layer.sum()) for layer in states]


def verify_exception_cases_moment3(
    cases: Set[Pair],
    case_metadata: Optional[Dict[int, Tuple[int, int]]] = None,
    quiet_dp: bool = False,
) -> bool:
    """
    Verify (q,k) cases by direct DP for (x,x^2,x^3), i.e. all q^3 triples.
    """
    if not cases:
        print("\nNo exceptional cases to verify with third-moment DP.")
        return True

    by_q: Dict[int, List[int]] = {}
    for q, k in cases:
        by_q.setdefault(q, []).append(k)

    print("\n" + "=" * 80)
    print("DIRECT DP VERIFICATION OF THIRD-MOMENT EXCEPTIONAL CASES")
    print("=" * 80)
    ok_all = True

    for q in sorted(by_q):
        ks = sorted(set(by_q[q]))
        if case_metadata and q in case_metadata:
            p, s = case_metadata[q]
        else:
            p, s = factor_prime_power(q)
        field = build_field_data(q, p, s)
        max_k = max(ks)
        target = q**3
        modulus_text = "prime field" if field.modulus is None else f"modulus={list(field.modulus)}"
        print(f"\nq={q}={p}^{s}, checking k={ks}, target size={target}, {modulus_text}")

        sizes = coverage_sizes_moment3_numpy(field, max_k, target_ks=set(ks), quiet=quiet_dp)
        for k in ks:
            ok = sizes[k] == target
            ok_all = ok_all and ok
            status = "OK" if ok else "FAIL"
            print(f"q={q:4d}, k={k:3d}: reached {sizes[k]:8d}/{target:8d} triples -- {status}")
    return ok_all


def format_ranges(ranges: Tuple[Tuple[int, int], ...]) -> str:
    if not ranges:
        return "-"
    return ", ".join(f"{{{a}}}" if a == b else f"[{a}, {b}]" for a, b in ranges)


def print_report(results: Sequence[CaseResult], failures: Set[Pair]) -> None:
    all_hold = [r for r in results if r.status == "all_hold"]
    some_fail = [r for r in results if r.status == "some_fail"]
    all_fail = [r for r in results if r.status == "all_fail"]
    empty = [r for r in results if r.status == "empty"]

    print("\n" + "=" * 80)
    print("RESULTS FOR THE NEW THIRD-MOMENT E(k,p,q) INEQUALITY")
    print("=" * 80)
    print(f"Total prime-power q: {len(results)}")
    print(f"All hold:            {len(all_hold)}")
    print(f"Some fail:           {len(some_fail)}")
    print(f"All fail:            {len(all_fail)}")
    print(f"Empty k-range:       {len(empty)}")

    print("\n" + "=" * 80)
    print(f"RUNTIME ANALYTIC FAILURE CASES  (total {len(failures)} points)")
    print("=" * 80)
    failures_by_case = sorted(some_fail + all_fail, key=lambda r: r.q)
    if not failures_by_case:
        print("None.")
    else:
        for r in failures_by_case:
            print(
                f"q={r.q:>5} = {r.p}^{r.s}, k in [{r.k_min}, {r.k_max}], "
                f"holds {r.hold_count}/{r.total_k}, fails at {format_ranges(r.fail_ranges)}"
            )


def compare_with_expected_analytic_failures(failures: Set[Pair]) -> bool:
    missing = EXPECTED_ANALYTIC_FAILURES - failures
    unexpected = failures - EXPECTED_ANALYTIC_FAILURES

    print("\n" + "=" * 80)
    print("REGRESSION CHECK AGAINST EXPECTED_ANALYTIC_FAILURES")
    print("=" * 80)
    print(f"Expected analytic failure points: {len(EXPECTED_ANALYTIC_FAILURES)}")
    print(f"Observed analytic failure points: {len(failures)}")

    if not missing and not unexpected:
        print("Observed analytic failures match EXPECTED_ANALYTIC_FAILURES.")
        return True

    if missing:
        print(f"Missing expected failures: {sorted(missing)}")
    if unexpected:
        print(f"Unexpected observed failures: {sorted(unexpected)}")
    return False


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Check the new third-moment E(k,p,q) inequality.")
    parser.add_argument("--qmin", type=int, default=Q_MIN)
    parser.add_argument("--qmax-exclusive", type=int, default=Q_MAX_EXCLUSIVE)
    parser.add_argument("--char-min", type=int, default=CHAR_MIN)
    parser.add_argument("--workers", type=int, default=cpu_count() or 1)
    parser.add_argument("--decimal", action="store_true", help="run selected high-precision Decimal checks")
    parser.add_argument(
        "--check-expected-analytic-failures",
        action="store_true",
        help="compare runtime analytic failures with EXPECTED_ANALYTIC_FAILURES",
    )
    parser.add_argument(
        "--verify-exceptions",
        action="store_true",
        help="directly verify EXPECTED_EXCEPTIONAL_CASES by third-moment DP",
    )
    parser.add_argument(
        "--verify-runtime-failures",
        action="store_true",
        help="directly verify every runtime analytic failure by third-moment DP",
    )
    parser.add_argument("--quiet-dp", action="store_true", help="suppress DP early-full progress lines")
    parser.add_argument("--verbose", action="store_true", help="print detailed data for each checked k")
    args = parser.parse_args(argv)

    cases = generate_prime_powers(args.qmin, args.qmax_exclusive, args.char_min)
    print("=" * 80)
    print("CHECKING NEW THIRD-MOMENT E(k,p,q) INEQUALITY")
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

    failures = {(r.q, k) for r in results for a, b in r.fail_ranges for k in range(a, b + 1)}
    case_metadata = {r.q: (r.p, r.s) for r in results}

    print_report(results, failures)

    if args.check_expected_analytic_failures:
        if not compare_with_expected_analytic_failures(failures):
            return 1

    if args.decimal:
        print("\n" + "=" * 80)
        print("SELECTED DECIMAL VERIFICATION")
        print("=" * 80)
        if verify_selected_points(results):
            return 1
        print("Selected Decimal verification passed.")

    if args.verify_exceptions:
        if not verify_exception_cases_moment3(EXPECTED_EXCEPTIONAL_CASES, case_metadata, quiet_dp=args.quiet_dp):
            print("ERROR: EXPECTED_EXCEPTIONAL_CASES did NOT pass third-moment DP verification!")
            return 1
        print("All EXPECTED_EXCEPTIONAL_CASES passed third-moment DP verification.")

    if args.verify_runtime_failures:
        if not verify_exception_cases_moment3(failures, case_metadata, quiet_dp=args.quiet_dp):
            print("ERROR: at least one runtime analytic failure does NOT pass third-moment DP verification.")
            return 1
        print("All runtime analytic failures passed third-moment DP verification.")

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
