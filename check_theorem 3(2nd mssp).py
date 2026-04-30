import math
import multiprocessing as mp
import sys
import argparse
from functools import lru_cache
from typing import List, Tuple

# ---------- Finite field (prime p) ----------
class GFp:
    """Integers modulo a prime p (p >= 5)."""
    def __init__(self, value: int, p: int):
        self.value = value % p
        self.p = p

    def __add__(self, other):
        return GFp(self.value + other.value, self.p)
    def __mul__(self, other):
        return GFp(self.value * other.value, self.p)
    def __sub__(self, other):
        return GFp(self.value - other.value, self.p)
    def __truediv__(self, other):
        inv = pow(other.value, self.p - 2, self.p)
        return GFp(self.value * inv, self.p)
    def __eq__(self, other):
        return self.value == other.value
    def __hash__(self):
        return hash((self.value, self.p))
    def __repr__(self):
        return str(self.value)
    def to_int(self) -> int:
        return self.value

# ---------- Combination unranking ----------
@lru_cache(maxsize=None)
def nCr(n: int, r: int) -> int:
    if r < 0 or r > n:
        return 0
    if r == 0 or r == n:
        return 1
    return nCr(n-1, r-1) + nCr(n-1, r)

def combination_from_rank(rank: int, n: int, k: int) -> List[int]:
    """Return the sorted list of k indices (0-based) for the given rank."""
    comb = []
    for i in range(k, 0, -1):
        x = i - 1
        while nCr(x, i) <= rank:
            x += 1
        comb.append(x - 1)
        rank -= nCr(x - 1, i)
    return comb

# ---------- Worker function ----------
def worker_task(args: Tuple):
    (p, k, start_rank, end_rank, q, worker_id) = args
    elements = [GFp(i, p) for i in range(p)]
    
    covered = [[False] * q for _ in range(q)]
    covered_count = 0
    target = q * q

    two = GFp(2, p)
    inv2 = GFp(1, p) / two

    n = p

    for rank in range(start_rank, end_rank):
        idx = combination_from_rank(rank, n, k)
        I = [elements[i] for i in idx]

        sum1 = GFp(0, p)
        sum_sq = GFp(0, p)
        for x in I:
            sum1 += x
            sum_sq += x * x
        sum2 = (sum1 * sum1 - sum_sq) * inv2

        u = sum1.to_int()
        v = sum2.to_int()
        if not covered[u][v]:
            covered[u][v] = True
            covered_count += 1
            if covered_count == target:
                return covered, True

    return covered, False

# ---------- Main verification for a single (q, k) ----------
def check_theorem(q: int, k: int, num_workers: int = None) -> bool:
    """Verify the theorem for GF(q), using multiple processes.
    Currently restricted to prime fields (q must be prime)."""
    if q < 17:
        raise ValueError("q must be >= 17")
    if not all(q % d for d in range(2, int(math.isqrt(q)) + 1)):
        raise NotImplementedError("Only prime fields are supported in this version. "
                                  "Use galois for extension fields.")

    if k < 4 or k > q - 4:
        print(f"q={q}, k={k}: k not in [4, q-4], skipping.")
        return True

    total_comb = nCr(q, k)
    print(f"q={q}, k={k}: total subsets = {total_comb}")

    if num_workers is None:
        num_workers = mp.cpu_count()
    num_workers = min(num_workers, total_comb)
    print(f"  Using {num_workers} workers")

    chunk_size = total_comb // num_workers
    ranges = []
    start = 0
    for i in range(num_workers):
        end = start + chunk_size + (1 if i < total_comb % num_workers else 0)
        ranges.append((start, end))
        start = end

    args = [(q, k, s, e, q, i) for i, (s, e) in enumerate(ranges)]

    with mp.Pool(processes=num_workers) as pool:
        results = pool.map(worker_task, args)

    covered = [[False] * q for _ in range(q)]
    all_done = False
    for local_covered, done in results:
        all_done = all_done or done
        for i in range(q):
            for j in range(q):
                if local_covered[i][j]:
                    covered[i][j] = True

    covered_count = sum(sum(1 for cell in row if cell) for row in covered)
    total_pairs = q * q

    if covered_count < total_pairs:
        missing = []
        for i in range(q):
            for j in range(q):
                if not covered[i][j]:
                    missing.append((GFp(i, q), GFp(j, q)))
                    if len(missing) >= 5:
                        break
            if len(missing) >= 5:
                break
        print(f"  FAILED: only {covered_count}/{total_pairs} pairs covered.")
        print(f"  First missing pairs (u,v): {missing}")
        return False
    else:
        print(f"  All {total_pairs} pairs covered.")
        return True

# ---------- CLI / Range mode ----------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Verify subset sum theorem over finite fields.")
    parser.add_argument("--k", type=int, required=True, help="Subset size")
    # Choice between a single q or a range of q:
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--q", type=int, default=None, help="Single prime field size")
    group.add_argument("--qmin", type=int, default=None, help="Minimum q for range test (with --qmax)")
    parser.add_argument("--qmax", type=int, default=None, help="Maximum q for range test")
    parser.add_argument("--workers", type=int, default=None, help="Number of processes")
    args = parser.parse_args()

    # Validate arguments
    if args.q is not None:
        # Single q mode
        q = args.q
        ok = check_theorem(q, args.k, args.workers)
        if ok:
            print("Verification successful.")
        else:
            print("Verification failed.")
    elif args.qmin is not None and args.qmax is not None:
        # Range mode
        if args.qmin > args.qmax:
            print("Error: qmin must be <= qmax")
            sys.exit(1)
        k = args.k
        print(f"Checking theorem for k={k} and prime q in [{args.qmin}, {args.qmax}] (char>=5, q>=17)")
        all_ok = True
        for q in range(max(args.qmin, 17), args.qmax + 1):
            # Only prime fields with characteristic >= 5
            if q < 5:
                continue
            if not all(q % d for d in range(2, int(math.isqrt(q)) + 1)):
                continue  # not prime
            if k < 4 or k > q - 4:
                print(f"q={q}: k out of range, skipping.")
                continue
            ok = check_theorem(q, k, args.workers)
            if not ok:
                all_ok = False
                print(f"Counterexample found at q={q}")
                break
        if all_ok:
            print("All tested q passed the theorem.")
    elif args.qmax is not None and args.qmin is None:
        # Legacy --qmax mode (without --qmin) : still works as before,
        # using qmin = 17 as default for backward compatibility
        k = args.k
        print(f"Checking theorem for k={k} and prime q up to {args.qmax} (starting from 17)")
        all_ok = True
        for q in range(17, args.qmax + 1):
            if not all(q % d for d in range(2, int(math.isqrt(q)) + 1)):
                continue
            if k < 4 or k > q - 4:
                print(f"q={q}: k out of range, skipping.")
                continue
            ok = check_theorem(q, k, args.workers)
            if not ok:
                all_ok = False
                print(f"Counterexample found at q={q}")
                break
        if all_ok:
            print("All tested q passed the theorem.")
    else:
        parser.error("You must specify either --q, or both --qmin and --qmax, or --qmax alone.")