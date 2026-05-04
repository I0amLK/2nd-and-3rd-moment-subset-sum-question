"""
Check the inequality for prime powers q with char(q) >= 5, 17 <= q <= 5381,
and integer k satisfying

    (1 + sqrt(8q - 7)) / 2 <= k <= q / 2.

The inequality is

    (q)_k / q^2
    - ((q^2 - q) / q^2)
      * (sqrt(q) + k - (q - sqrt(q))/p - 1)_k
    > 0.

Equivalently, since q^2 > 0,

    (q)_k > (q^2 - q) * (a_k)_k,

where

    a_k = sqrt(q) + k - (q - sqrt(q))/p - 1.

We use the recurrence

    f(k) = (q)_k / ((q^2 - q) * (a_k)_k)

and

    f(k+1) / f(k)
    = (q - k) / (sqrt(q) + k - (q - sqrt(q))/p).

Therefore

    log|f(k+1)|
    =
    log|f(k)| + log(q - k) - log|d_k|,

where

    d_k = sqrt(q) + k - (q - sqrt(q))/p.

The sign is updated by multiplying by sign(d_k).
"""

import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from decimal import Decimal, getcontext
from os import cpu_count
from typing import Dict, List, Optional, Tuple

from tqdm import tqdm


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

Q_MIN = 17
Q_MAX = 5381
CHAR_MIN = 5

ZERO_TOL = 1e-14


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Prime and prime-power utilities
# ---------------------------------------------------------------------------

def is_prime(n: int) -> bool:
    """Return True iff n is prime."""
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


def generate_prime_powers(
    lo: int = Q_MIN,
    hi: int = Q_MAX,
    char_min: int = CHAR_MIN,
) -> List[Tuple[int, int, int]]:
    """
    Return sorted list of (q, p, s), where q = p^s,
    p is prime, and p >= char_min.
    """
    cases = []
    seen = set()

    for p in range(char_min, hi + 1):
        if not is_prime(p):
            continue

        q = p
        s = 1

        while q <= hi:
            if q >= lo and q not in seen:
                seen.add(q)
                cases.append((q, p, s))

            q *= p
            s += 1

    cases.sort()
    return cases


# ---------------------------------------------------------------------------
# k-range utilities
# ---------------------------------------------------------------------------

def k_lower(q: int) -> int:
    return math.ceil((1.0 + math.sqrt(8.0 * q - 7.0)) / 2.0)


def k_upper(q: int) -> int:
    return q // 2


# ---------------------------------------------------------------------------
# Failure-range utility
# ---------------------------------------------------------------------------

def collapse_to_ranges(values: List[int]) -> Tuple[Tuple[int, int], ...]:
    """Collapse [1,2,3,7,8] to ((1,3), (7,8))."""
    if not values:
        return tuple()

    values = sorted(values)
    ranges = []

    start = values[0]
    end = values[0]

    for x in values[1:]:
        if x == end + 1:
            end = x
        else:
            ranges.append((start, end))
            start = end = x

    ranges.append((start, end))
    return tuple(ranges)


# ---------------------------------------------------------------------------
# Core recurrence check
# ---------------------------------------------------------------------------

def check_inequality_range(
    q: int,
    p: int,
    s: int,
    verbose: bool = False,
) -> Dict[int, bool]:
    """
    Return a dictionary mapping k -> whether the inequality holds.

    The check is done by computing f(k_min) once, and then updating
    f(k) recursively.
    """
    sqrt_q = math.sqrt(q)

    k_min = k_lower(q)
    k_max = k_upper(q)

    if k_min > k_max:
        return {}

    # Let
    #
    #   offset = sqrt(q) - (q - sqrt(q)) / p.
    #
    # Then
    #
    #   a_k = sqrt(q) + k - (q - sqrt(q))/p - 1
    #       = k + offset - 1.
    #
    offset = sqrt_q - (q - sqrt_q) / p

    # ------------------------------------------------------------
    # Initialize g = log|f(k_min)|.
    #
    # f(k) = (q)_k / ((q^2 - q) * (a_k)_k).
    # ------------------------------------------------------------

    ln_q_ff = 0.0
    for i in range(k_min):
        ln_q_ff += math.log(q - i)

    a_init = k_min + offset - 1.0

    ln_a_ff = 0.0

    # Since (q)_k > 0 and q^2 - q > 0,
    # sign(f(k)) equals sign((a_k)_k).
    sign_f = 1

    for i in range(k_min):
        val = a_init - i

        if abs(val) < ZERO_TOL:
            # Then (a_k)_k = 0, so f is formally infinite.
            ln_a_ff = float("-inf")
            sign_f = 0
            break

        if val < 0:
            sign_f = -sign_f

        ln_a_ff += math.log(abs(val))

    ln_q2_minus_q = math.log(q * (q - 1))

    if sign_f == 0:
        g = float("inf")
    else:
        g = ln_q_ff - ln_q2_minus_q - ln_a_ff

    if verbose:
        print(
            f"Initial q={q}, p={p}, s={s}, "
            f"k_min={k_min}, k_max={k_max}, "
            f"offset={offset:.12g}, g={g:.12g}, sign_f={sign_f}"
        )

    # ------------------------------------------------------------
    # Iterate through all k.
    # ------------------------------------------------------------

    results: Dict[int, bool] = {}

    for k in range(k_min, k_max + 1):
        # Check the sign and magnitude.
        #
        # If (a_k)_k <= 0, then RHS <= 0, while LHS = (q)_k > 0,
        # so the inequality automatically holds.
        #
        # If (a_k)_k > 0, we need f(k) > 1, i.e. log(f(k)) > 0.
        if sign_f <= 0:
            holds = True
        else:
            holds = g > 0.0

        results[k] = holds

        if verbose and k <= k_min + 5:
            print(
                f"  k={k}, g={g:.12g}, sign_f={sign_f}, holds={holds}"
            )

        # Update from k to k + 1.
        if k < k_max:
            d_k = k + offset
            numerator_factor = q - k

            if numerator_factor <= 0:
                raise ValueError(f"Unexpected q-k <= 0 for q={q}, k={k}")

            if abs(d_k) < ZERO_TOL:
                # The denominator falling factorial at the next k has
                # a zero factor, hence RHS is zero and inequality holds
                # at that next k. We mark sign_f = 0.
                g = float("inf")
                sign_f = 0
            else:
                g = g + math.log(numerator_factor) - math.log(abs(d_k))

                if sign_f == 0:
                    # In exact arithmetic, after a zero factor appears,
                    # the next product may no longer be recoverable by this
                    # ratio because the previous product was zero.
                    #
                    # This case is extremely rare and only happens when d_k
                    # is exactly zero in floating approximation. To be safe,
                    # recompute the next sign and log directly.
                    next_k = k + 1
                    sign_f, g = recompute_sign_and_log_f(q, p, next_k)
                elif d_k < 0:
                    sign_f = -sign_f

    return results


def recompute_sign_and_log_f(q: int, p: int, k: int) -> Tuple[int, float]:
    """
    Recompute sign(f(k)) and log|f(k)| directly.

    Used only as a safety fallback near zero-factor cases.
    """
    sqrt_q = math.sqrt(q)
    offset = sqrt_q - (q - sqrt_q) / p
    a_k = k + offset - 1.0

    ln_q_ff = 0.0
    for i in range(k):
        ln_q_ff += math.log(q - i)

    ln_a_ff = 0.0
    sign_f = 1

    for i in range(k):
        val = a_k - i

        if abs(val) < ZERO_TOL:
            return 0, float("inf")

        if val < 0:
            sign_f = -sign_f

        ln_a_ff += math.log(abs(val))

    g = ln_q_ff - math.log(q * (q - 1)) - ln_a_ff
    return sign_f, g


# ---------------------------------------------------------------------------
# Per-case worker for multiprocessing
# ---------------------------------------------------------------------------

def check_one_case(case: Tuple[int, int, int]) -> CaseResult:
    q, p, s = case

    k_min = k_lower(q)
    k_max = k_upper(q)

    if k_min > k_max:
        return CaseResult(
            q=q,
            p=p,
            s=s,
            status="empty",
            k_min=None,
            k_max=None,
            total_k=0,
            hold_count=0,
            fail_count=0,
            fail_ranges=tuple(),
        )

    results = check_inequality_range(q, p, s)

    total_k = len(results)
    hold_count = sum(1 for v in results.values() if v)
    fail_count = total_k - hold_count

    fail_ks = [k for k, ok in results.items() if not ok]
    fail_ranges = collapse_to_ranges(fail_ks)

    if fail_count == 0:
        status = "all_hold"
    elif hold_count == 0:
        status = "all_fail"
    else:
        status = "some_fail"

    return CaseResult(
        q=q,
        p=p,
        s=s,
        status=status,
        k_min=k_min,
        k_max=k_max,
        total_k=total_k,
        hold_count=hold_count,
        fail_count=fail_count,
        fail_ranges=fail_ranges,
    )


# ---------------------------------------------------------------------------
# Optional Decimal verification
# ---------------------------------------------------------------------------

def decimal_bruteforce_check(q: int, p: int, k: int, prec: int = 200) -> bool:
    """
    Brute-force Decimal check for one q,p,k.

    This is slower, so use it only for verification.
    """
    getcontext().prec = prec

    qd = Decimal(q)
    pd = Decimal(p)
    kd = Decimal(k)

    sqrt_qd = qd.sqrt()

    inner = sqrt_qd + kd - (qd - sqrt_qd) / pd - Decimal(1)

    def falling_decimal(x: Decimal, n: int) -> Decimal:
        r = Decimal(1)
        for i in range(n):
            r *= x - Decimal(i)
        return r

    lhs = falling_decimal(qd, k)
    rhs = (qd * qd - qd) * falling_decimal(inner, k)

    return lhs > rhs


def verify_selected_points(results: List[CaseResult]) -> int:
    """
    Verify selected k values with Decimal arithmetic.

    For each non-empty case, verify k_min, midpoint, and k_max.
    For cases with failures, also verify endpoints of failure ranges.
    """
    mismatches = 0
    checked = 0

    for r in tqdm(results, desc="Decimal verification", unit="case"):
        if r.status == "empty":
            continue

        assert r.k_min is not None
        assert r.k_max is not None

        recurrence_results = check_inequality_range(r.q, r.p, r.s)

        ks = {r.k_min, r.k_max, (r.k_min + r.k_max) // 2}

        for a, b in r.fail_ranges:
            ks.add(a)
            ks.add(b)
            if a > r.k_min:
                ks.add(a - 1)
            if b < r.k_max:
                ks.add(b + 1)

        for k in sorted(ks):
            expected = recurrence_results[k]
            actual = decimal_bruteforce_check(r.q, r.p, k)
            checked += 1

            if expected != actual:
                print(
                    f"MISMATCH: q={r.q}={r.p}^{r.s}, k={k}, "
                    f"recurrence={expected}, decimal={actual}"
                )
                mismatches += 1

    print(f"Decimal verification checked {checked} selected k-values.")
    return mismatches


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def format_ranges(ranges: Tuple[Tuple[int, int], ...]) -> str:
    if not ranges:
        return "-"

    parts = []
    for a, b in ranges:
        if a == b:
            parts.append(f"{{{a}}}")
        else:
            parts.append(f"[{a}, {b}]")

    return ", ".join(parts)


def print_report(results: List[CaseResult]) -> None:
    all_hold = [r for r in results if r.status == "all_hold"]
    some_fail = [r for r in results if r.status == "some_fail"]
    all_fail = [r for r in results if r.status == "all_fail"]
    empty = [r for r in results if r.status == "empty"]

    print()
    print("=" * 80)
    print("RESULTS")
    print("=" * 80)

    print(f"All hold:   {len(all_hold)}")
    print(f"Some fail:  {len(some_fail)}")
    print(f"All fail:   {len(all_fail)}")
    print(f"Empty:      {len(empty)}")

    print()
    print("=" * 80)
    print("SOME FAIL CASES")
    print("=" * 80)

    if not some_fail:
        print("None.")
    else:
        for r in some_fail:
            print(
                f"q={r.q:>4} = {r.p}^{r.s}, "
                f"k in [{r.k_min}, {r.k_max}], "
                f"holds {r.hold_count}/{r.total_k}, "
                f"fails at k in {format_ranges(r.fail_ranges)}"
            )

    print()
    print("=" * 80)
    print("ALL FAIL CASES")
    print("=" * 80)

    if not all_fail:
        print("None.")
    else:
        for r in all_fail:
            print(
                f"q={r.q:>4} = {r.p}^{r.s}, "
                f"k in [{r.k_min}, {r.k_max}], "
                f"total k={r.total_k}"
            )

    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)

    total = len(results)

    def pct(n: int) -> float:
        return 100.0 * n / total if total else 0.0

    print(f"Total cases: {total}")
    print(f"All hold:    {len(all_hold)} ({pct(len(all_hold)):.2f}%)")
    print(f"Some fail:   {len(some_fail)} ({pct(len(some_fail)):.2f}%)")
    print(f"All fail:    {len(all_fail)} ({pct(len(all_fail)):.2f}%)")
    print(f"Empty:       {len(empty)} ({pct(len(empty)):.2f}%)")

    primes = [r for r in results if r.s == 1]
    proper_pp = [r for r in results if r.s >= 2]

    print()
    print("By exponent:")
    print(f"Primes s=1:               {len(primes)}")
    print(f"Proper prime powers s>=2: {len(proper_pp)}")

    print()
    print("Prime cases:")
    print(f"  all hold:  {sum(1 for r in primes if r.status == 'all_hold')}")
    print(f"  some fail: {sum(1 for r in primes if r.status == 'some_fail')}")
    print(f"  all fail:  {sum(1 for r in primes if r.status == 'all_fail')}")
    print(f"  empty:     {sum(1 for r in primes if r.status == 'empty')}")

    print()
    print("Proper prime-power cases:")
    print(f"  all hold:  {sum(1 for r in proper_pp if r.status == 'all_hold')}")
    print(f"  some fail: {sum(1 for r in proper_pp if r.status == 'some_fail')}")
    print(f"  all fail:  {sum(1 for r in proper_pp if r.status == 'all_fail')}")
    print(f"  empty:     {sum(1 for r in proper_pp if r.status == 'empty')}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=" * 80)
    print("CHECKING INEQUALITY")
    print("=" * 80)
    print(f"q range:      [{Q_MIN}, {Q_MAX}]")
    print(f"char(q) >=    {CHAR_MIN}")

    cases = generate_prime_powers(Q_MIN, Q_MAX, CHAR_MIN)

    print(f"Prime-power cases found: {len(cases)}")

    workers = cpu_count() or 4
    print(f"Workers: {workers}")

    results: List[CaseResult] = []

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(check_one_case, case) for case in cases]

        for future in tqdm(as_completed(futures), total=len(futures), desc="Checking q", unit="case"):
            try:
                results.append(future.result())
            except Exception as exc:
                print(f"Worker error: {exc}")
                raise

    results.sort(key=lambda r: r.q)

    print_report(results)

    # Optional verification.
    #
    # This is slower than the recurrence. Leave it on if you want strong
    # confidence; turn it off if you only want the fast run.
    DO_DECIMAL_VERIFICATION = True

    if DO_DECIMAL_VERIFICATION:
        print()
        print("=" * 80)
        print("DECIMAL VERIFICATION")
        print("=" * 80)

        mismatches = verify_selected_points(results)

        if mismatches:
            print(f"Decimal verification found {mismatches} mismatches.")
            sys.exit(1)

        print("Decimal verification passed.")

    print()
    print("Done.")


if __name__ == "__main__":
    main()