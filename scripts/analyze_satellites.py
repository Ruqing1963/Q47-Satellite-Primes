#!/usr/bin/env python3
"""
analyze_satellites.py — Reproduce all statistical analyses from Paper IV
(Satellite Primes: The Local Prime Landscape around Giant Primes).

Reads the satellite catalog and computes:
  1. Poisson fit with zero-satellite recovery (N=2107)
  2. Gap uniformity chi-square test
  3. Mod-30 residue structure
  4. Nearest-satellite CDF vs Cramér exponential
  5. Conditional Hardy-Littlewood analysis (Bayesian doubling)
  6. Density vs n comparison

Usage:
    python analyze_satellites.py [--data ../data/satellites_9012.csv]

Author: Ruqing Chen, GUT Geoservice Inc.
Date:   February 2026
"""

import argparse
import csv
import numpy as np
from collections import Counter, defaultdict
from math import factorial, exp
from scipy import stats


def load_satellites(path):
    """Load satellite catalog from CSV."""
    satellites = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            satellites.append((int(row['main_star_n']), int(row['gap_k'])))
    return satellites


def compute_S(k, max_p=5000):
    """Compute unconditional HL singular series for pair (n, n+k)."""
    product = 1.0
    sieve = [True] * (max_p + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(max_p**0.5) + 1):
        if sieve[i]:
            for j in range(i * i, max_p + 1, i):
                sieve[j] = False
    for p in range(2, max_p + 1):
        if not sieve[p]:
            continue
        bad = set()
        for r in range(p):
            if r % p == 0:
                bad.add(r)
            if (r + k) % p == 0:
                bad.add(r)
        nu = len(bad)
        product *= (1 - nu / p) / (1 - 1 / p) ** 2
    return product


def main():
    parser = argparse.ArgumentParser(
        description='Reproduce Paper IV satellite analyses')
    parser.add_argument('--data', '-d', default='../data/satellites_9012.csv',
                        help='Path to satellite catalog CSV')
    args = parser.parse_args()

    satellites = load_satellites(args.data)
    by_star = defaultdict(list)
    for n, k in satellites:
        by_star[n].append(k)
    stars = sorted(by_star.keys())
    all_gaps = [k for _, k in satellites]
    sats_per = [len(by_star[n]) for n in stars]
    n_with = len(stars)

    print("=" * 70)
    print("  PAPER IV: SATELLITE PRIME ANALYSIS")
    print("=" * 70)

    # 1. Poisson fit with zero-satellite recovery
    print("\n--- 1. POISSON FIT (CORRECTED) ---")
    lam_naive = np.mean(sats_per)
    p_zero = exp(-lam_naive)
    N_true = int(round(n_with / (1 - p_zero)))
    N_zero = N_true - n_with
    lam_corrected = len(satellites) / N_true

    print(f"  Stars with ≥1 satellite: {n_with}")
    print(f"  Inferred total (N_true): {N_true}")
    print(f"  Zero-satellite stars: ~{N_zero}")
    print(f"  Corrected λ = {len(satellites)}/{N_true} = {lam_corrected:.3f}")

    all_sats = [0] * N_zero + sats_per
    disp = np.var(all_sats) / np.mean(all_sats)
    print(f"  Dispersion index: {disp:.3f}")

    print(f"\n  {'k':>3} {'Obs':>6} {'Poisson':>8} {'Ratio':>6}")
    for k in range(15):
        obs = N_zero if k == 0 else sum(1 for s in sats_per if s == k)
        poi = N_true * (lam_corrected ** k * exp(-lam_corrected) / factorial(k))
        print(f"  {k:>3} {obs:>6} {poi:>8.1f} {obs / poi:>6.2f}")

    # 2. Gap uniformity
    print("\n--- 2. GAP UNIFORMITY ---")
    bins = np.arange(0, 5500, 500)
    hist, _ = np.histogram(all_gaps, bins=bins)
    chi2, pval = stats.chisquare(hist)
    print(f"  Chi-square (10 bins): {chi2:.2f}, p-value: {pval:.4f}")

    # 3. Mod-30 structure
    print("\n--- 3. MOD-30 RESIDUE STRUCTURE ---")
    mod30 = Counter(k % 30 for k in all_gaps)
    admissible = [0, 2, 6, 8, 12, 14, 18, 20, 24, 26]
    print(f"  {'k%30':>5} {'k%6':>4} {'Count':>6} {'%':>6}")
    for r in admissible:
        print(f"  {r:>5} {r % 6:>4} {mod30[r]:>6} {mod30[r] / len(all_gaps) * 100:>5.1f}%")

    # 4. Nearest-satellite CDF
    print("\n--- 4. NEAREST-SATELLITE CDF ---")
    min_gaps = [min(by_star[n]) for n in stars]
    ln_P = np.mean([int(46 * np.log10(n) + 1.67) * np.log(10) for n in stars])
    print(f"  Average ln(P): {ln_P:.1f}")
    for thresh in [50, 100, 200, 500, 1000, 2000, 3000]:
        obs_cdf = sum(1 for g in min_gaps if g <= thresh) / n_with
        cramer_cdf = 1 - exp(-thresh / (3 * ln_P))
        print(f"  k ≤ {thresh:>4}: Cramér={cramer_cdf:.3f}, Obs={obs_cdf:.3f}, "
              f"ratio={obs_cdf / cramer_cdf:.2f}")

    # 5. Conditional HL (Bayesian doubling)
    print("\n--- 5. CONDITIONAL HARDY-LITTLEWOOD ---")
    gap_counts = Counter(all_gaps)
    print(f"  {'k':>3} {'k%6':>4} {'B(k)':>5} {'S_cond':>7} {'E':>6} {'Obs':>4}")
    for k in [2, 6, 8, 12, 14, 18, 20, 24, 26, 30]:
        S = compute_S(k)
        B = 2 if k % 3 == 2 else 1
        S_cond = S * B
        E = N_true * S_cond / ln_P
        obs = gap_counts.get(k, 0)
        print(f"  {k:>3} {k % 6:>4} {B:>5} {S_cond:>7.3f} {E:>6.2f} {obs:>4}")

    # 6. Density vs n
    print("\n--- 6. SATELLITE DENSITY VS n ---")
    n_arr = np.array(stars)
    s_arr = np.array(sats_per)
    for lo, hi in [(50e9, 75e9), (75e9, 100e9), (100e9, 125e9),
                   (125e9, 150e9), (150e9, 175e9), (175e9, 200e9)]:
        mask = (n_arr >= lo) & (n_arr < hi)
        if mask.sum() == 0:
            continue
        d_mean = np.mean(46 * np.log10(n_arr[mask]) + 1.67)
        cramer = 5000 / (d_mean * np.log(10))
        obs_mean = np.mean(s_arr[mask])
        print(f"  [{lo / 1e9:.0f}B,{hi / 1e9:.0f}B): {mask.sum():>4} stars, "
              f"obs={obs_mean:.2f}, Cramér={cramer:.2f}, "
              f"ratio={obs_mean / cramer:.3f}")

    print("\n" + "=" * 70)
    print("  ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
