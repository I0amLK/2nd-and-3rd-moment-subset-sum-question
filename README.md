# 2nd and 3rd Moment Subset-Sum Question

This repository contains the programs used to study the **second- and third-moment
subset-sum problems** over finite fields.

The core question is: *for which `q` and `k` does the family of all `k`-subsets
of `GF(q)` attain every possible value of the moment tuple*

```
(Σ x,  Σ x²,  Σ x³)   in   GF(q)³
```

(and the second-moment variant, the pair `(Σ x, Σ x²)`)?

The scripts fall into two families:

* **Brute-force / DP coverage checks** — directly enumerate `k`-subsets (or run a
  bitmask dynamic program) and verify whether every moment tuple is attained.
* **Analytic inequality checks** — verify the large-range lower-bound inequality
  `E(k, p, q) > 0` that is expected to hold away from the small exceptional cases.

Each script is self-contained and has no third-party dependencies beyond a
standard Python install (the `.py` files) or a Magma installation (the `.m`
files).

---

## File overview

| File | Language | Moment | Method | Description |
|------|----------|--------|--------|-------------|
| `check_2nd_moment_single.m` | Magma | 2nd | Brute force | Check one `(q, k)` pair. |
| `check_2nd_moment_range.m` | Magma | 2nd | Brute force | Check a fixed `k` over a range of `q`. |
| `check_2nd_moment_range.py` | Python | 2nd | Brute force (multiprocess) | Check a fixed `k` over a range of `q`. |
| `check_3rd_moment_small_q.m` | Magma | 3rd | Bitmask DP | Fast check of all triples for small `q`. |
| `check_2nd_moment_inequality.py` | Python | 2nd | Analytic inequality | Verify `E(k, q, p) > 0` over a large range. |
| `check_3rd_moment_inequality.py` | Python | 3rd | Analytic inequality + DP | Verify the new `E(k, p, q) > 0`, with finite-field DP for exceptional cases. |

---

## File descriptions

### `check_2nd_moment_single.m`
The simplest checker. Defines `CheckTheorem(q, k)`, which enumerates all
`k`-subsets of `GF(q)` and tests whether every pair `(Σ x, Σ x²)` is attained.
Returns `true` if all `q²` pairs are covered, otherwise `false` and prints the
first missing pairs. The file ends with the example call `CheckTheorem(17, 4)`.

> **Note:** requires `char(GF(q)) ≠ 2` because it divides by 2 when computing
> `Σ_{i<j} x_i x_j = ((Σ x)² − Σ x²) / 2`.

### `check_2nd_moment_range.m`
Extended version of the single-pair checker. It loops `q` from `17` up to a
configurable `Qmax` (prime powers of characteristic `≥ 5` only) for a fixed
subset size `k` (default `6`), applying the admissibility range `4 ≤ k ≤ q − 4`.
Stops on the first failing `q`. Edit the parameters block (`k`, `Qmax`) at the
top of the file to configure a run.

### `check_2nd_moment_range.py`
Python port of the range checker above, parallelised with `multiprocessing`.
For a fixed `k` it splits the `binom(q, k)` subsets into rank ranges handled by
worker processes, and reports — for each prime `q` in the requested range —
whether all pairs `(Σ x, Σ x²)` are covered. This is the brute-force second-
moment subset-sum (2nd MSSP) checker.

### `check_3rd_moment_small_q.m`
Fast third-moment checker for small fields. Instead of enumerating subsets, it
builds a bitmask dynamic program: `DP[t+1][(s1,s2)]` is a `q`-bit mask whose
`s3` bit is set iff some `t`-subset has moments `(s1, s2, s3)`. It then checks
whether all triples `(Σ x, Σ x², Σ x³)` are attained. Supports prime fields and
degree-2 extensions, uses complement symmetry (`k ≤ q/2`) and a centered
`s1 = 0` reduction when `p ∤ k`. Default range: `17 ≤ q ≤ 50`, `char ≥ 5`,
`k ≥ 6`.

### `check_2nd_moment_inequality.py`
Analytic checker for the corrected large-range second-moment inequality

```
E(k, q, p) = (q)_k − (q − 1)·L_{k,q,p} − (q² − q)·Q_{k,q,p} > 0,
```

where `(q)_k` is the falling factorial, `L_{k,q,p} = k!·C(q/p + k/p − 1, k/p)`
when `p | k` (else `0`), and `Q_{k,q,p} = (√q + k + (q − √q)/p − 1)_k`. It scans
prime powers `q = p^s` with `p ≥ 5` and `17 ≤ q < 53267` over the admissible
range `(1 + √(8q − 7))/2 < k ≤ q/2`, using multiprocessing and optional
high-precision `Decimal` arithmetic.

### `check_3rd_moment_inequality.py`
Analytic checker for the **new** third-moment lower-bound inequality

```
E(k, p, q) = (q)_k/q³ − ((q³ − q)/q³)·(L_{k,q,p} + Q_{k,q,p}) > 0,
```

equivalently `(q)_k > (q³ − q)·(L_{k,q,p} + Q_{k,q,p})`, with `L` and `Q` as
above. For the small *expected exceptional cases* where the analytic bound is
not expected to hold, it falls back to a direct finite-field DP check of the
third-moment subset-sum question (coverage of all triples
`(Σ x, Σ x², Σ x³)` over `GF(p^s)`, using NumPy boolean arrays). Default range
matches the second-moment inequality script.

---

## How to use

### Prerequisites
* **Python 3.8+** for the `.py` scripts. `check_3rd_moment_inequality.py` also
  uses NumPy for its DP fallback. Install with `pip install numpy` if needed.
* **Magma** for the `.m` scripts.

### Magma scripts

Run any `.m` file directly:

```bash
magma check_2nd_moment_single.m
magma check_2nd_moment_range.m
magma check_3rd_moment_small_q.m
```

`check_2nd_moment_single.m` defines a function — to check a different case, edit
the example call at the bottom of the file:

```magma
CheckTheorem(17, 4);   // change (q, k) here
```

`check_2nd_moment_range.m` and `check_3rd_moment_small_q.m` expose a parameters
block at the top of the file (e.g. `k`, `Qmax`, `MaxQ`, `MinK`, `UseSymmetry`).
Edit those values to configure a run, then launch with `magma <file>`.

### Python scripts

**`check_2nd_moment_range.py`** — brute-force 2nd-moment check for a fixed `k`
over a range of prime `q`:

```bash
python check_2nd_moment_range.py --k 4 --qmin 23 --qmax 43 --workers 8
```

Check a single `q` instead of a range:

```bash
python check_2nd_moment_range.py --k 4 --q 29
```

| Flag | Description |
|------|-------------|
| `--k` | Subset size (required). |
| `--q` | Test a single prime field size. |
| `--qmin` | Minimum `q` for a range run (use with `--qmax`). |
| `--qmax` | Maximum `q` for a range run. |
| `--workers` | Number of worker processes (defaults to all CPUs). |

**`check_2nd_moment_inequality.py`** — analytic 2nd-moment inequality:

```bash
python check_2nd_moment_inequality.py
python check_2nd_moment_inequality.py --qmin 17 --qmax-exclusive 53267 --workers 8 --verbose
python check_2nd_moment_inequality.py --decimal --verify-exceptions
```

**`check_3rd_moment_inequality.py`** — analytic 3rd-moment inequality with DP
fallback for exceptional cases:

```bash
python check_3rd_moment_inequality.py
python check_3rd_moment_inequality.py --qmin 17 --qmax-exclusive 53267 --workers 8 --verbose
python check_3rd_moment_inequality.py --decimal --quiet-dp
```

Common flags for the two inequality scripts:

| Flag | Description |
|------|-------------|
| `--qmin` | Smallest `q` to check (default `17`). |
| `--qmax-exclusive` | Exclusive upper bound on `q` (default `53267`). |
| `--char-min` | Smallest characteristic `p` (default `5`). |
| `--workers` | Number of worker processes (default: all CPUs). |
| `--decimal` | Run selected checks with high-precision `Decimal` arithmetic. |
| `--verbose` | Print detailed data for each checked `k`. |

Additional flags for `check_3rd_moment_inequality.py`: `--verify-exceptions`
forces the finite-field DP check on the expected exceptional cases, and
`--quiet-dp` suppresses DP progress lines.

---

## Typical workflow

1. **Quick sanity check** of a single case with
   `magma check_2nd_moment_single.m` (edit the `CheckTheorem(q, k)` call).
2. **Brute-force scan** of small `q` for a fixed `k` with
   `check_2nd_moment_range.py` (Python, parallel) or
   `check_2nd_moment_range.m` (Magma).
3. **Third-moment coverage** for small `q` with
   `magma check_3rd_moment_small_q.m` (fast DP).
4. **Large-range verification** of the analytic bounds with
   `check_2nd_moment_inequality.py` and `check_3rd_moment_inequality.py`, which
   cover `q` up to `53267` and fall back to a direct DP check on the small
   exceptional cases that the analytic inequality does not cover.
