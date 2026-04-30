# 2nd and 3rd Moment Subset-Sum Question
This repository contains programs used to study the second- and third-moment subset-sum problems.

## Usage

### check_theorem.m
Call `CheckTheorem(q, k)` to verify whether the theorem holds for the given `q` and `k`.

### check_theorem 3 (2nd MSSP).py
This Python script tests the theorem for a fixed `k` across a range of `q` values.

Example

Run:
`python theorem_check.py --k 4 --qmin 23 --qmax 43 --workers 8`

This command checks all cases with `k = 4` and `q` in `[23, 43]`.