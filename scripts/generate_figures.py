#!/usr/bin/env python3
"""
generate_figures.py — Reproduce all 4 figures for Paper IV
(Satellite Primes: The Local Prime Landscape).

Requires: matplotlib, numpy, scipy
Data:     ../data/satellites_9012.csv

Usage:
    python generate_figures.py

Author: Ruqing Chen, GUT Geoservice Inc.
Date:   February 2026
"""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
from math import factorial, exp

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 11, 'axes.labelsize': 12,
    'axes.titlesize': 13, 'figure.dpi': 300, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'mathtext.fontset': 'cm',
})

DATA = os.path.join(os.path.dirname(__file__), '..', 'data')
OUT = os.path.join(os.path.dirname(__file__), '..', 'paper', 'figures')
os.makedirs(OUT, exist_ok=True)

# Load data
satellites = []
with open(os.path.join(DATA, 'satellites_9012.csv')) as f:
    reader = csv.DictReader(f)
    for row in reader:
        satellites.append((int(row['main_star_n']), int(row['gap_k'])))

by_star = defaultdict(list)
for n, k in satellites:
    by_star[n].append(k)
stars = sorted(by_star.keys())
all_gaps = [k for _, k in satellites]
sats_per = [len(by_star[n]) for n in stars]
n_with = len(stars)

# Corrected parameters
N_true = 2107
N_zero = N_true - n_with
lam = 9012 / N_true  # 4.278


# ── Figure 1: Gap distribution + mod 30 ──
def figure1():
    from matplotlib.patches import Patch
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    bins_hist = np.arange(0, 5250, 250)
    ax1.hist(all_gaps, bins=bins_hist, color='#3498db', alpha=0.7,
             edgecolor='#2c3e50', linewidth=0.5, label='Observed')
    uniform = len(all_gaps) * 250 / 5000
    ax1.axhline(y=uniform, color='#e74c3c', linestyle='--', linewidth=1.5,
                label=f'Uniform ({uniform:.0f}/bin)')
    ax1.set_xlabel(r'Gap $k$ in $P - k$'); ax1.set_ylabel('Count')
    ax1.set_title('(a) Gap distribution (uniform, $p = 0.31$)', fontweight='bold')
    ax1.legend(fontsize=10); ax1.grid(True, alpha=0.15, axis='y')

    mod30 = Counter(k % 30 for k in all_gaps)
    admissible = [0, 2, 6, 8, 12, 14, 18, 20, 24, 26]
    vals = [mod30.get(r, 0) for r in admissible]
    colors = ['#e74c3c' if r % 6 == 0 else '#3498db' for r in admissible]
    ax2.bar(range(10), vals, color=colors, alpha=0.8, edgecolor='#333', linewidth=0.5)
    ax2.set_xticks(range(10))
    ax2.set_xticklabels([str(r) for r in admissible], fontsize=9)
    ax2.set_xlabel('$k$ mod 30'); ax2.set_ylabel('Count')
    ax2.set_title('(b) Mod-30 residue structure', fontweight='bold')
    ax2.axhline(y=len(all_gaps) / 10, color='grey', linestyle=':', alpha=0.5)
    ax2.legend([Patch(facecolor='#e74c3c', alpha=0.8),
                Patch(facecolor='#3498db', alpha=0.8)],
               [r'$k \equiv 0$ (mod 6)', r'$k \equiv 2$ (mod 6)'], fontsize=9)
    ax2.grid(True, alpha=0.15, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, 'p3_fig1.png'), dpi=300); plt.close()
    print("  Figure 1 saved.")


# ── Figure 2: CDF + Poisson ──
def figure2():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    min_gaps = sorted([min(by_star[n]) for n in stars])
    cdf_obs = np.arange(1, n_with + 1) / n_with
    ln_P = np.mean([int(46 * np.log10(n) + 1.67) * np.log(10) for n in stars])
    k_th = np.linspace(2, 5000, 500)
    cdf_cr = 1 - np.exp(-k_th / (3 * ln_P))
    ax1.plot(min_gaps, cdf_obs, '-', color='#3498db', linewidth=2, label='Observed CDF')
    ax1.plot(k_th, cdf_cr, '--', color='#e74c3c', linewidth=1.5,
             label=r"Cram$\acute{\mathrm{e}}$r: $1 - e^{-k/3\ln P}$")
    ax1.set_xlabel('Nearest satellite gap $k$'); ax1.set_ylabel('Cumulative probability')
    ax1.set_title('(a) Nearest satellite CDF', fontweight='bold')
    ax1.legend(fontsize=10); ax1.grid(True, alpha=0.2); ax1.set_xlim(0, 5000)

    k_range = range(16)
    obs_hist = [N_zero if k == 0 else sum(1 for s in sats_per if s == k) for k in k_range]
    poi_exp = [N_true * (lam ** k * exp(-lam) / factorial(k)) for k in k_range]
    x = np.arange(len(k_range)); w = 0.35
    ax2.bar(x - w / 2, obs_hist, w, color='#3498db', alpha=0.8, label='Observed')
    ax2.bar(x + w / 2, poi_exp, w, color='#e74c3c', alpha=0.6,
            label=f'Poisson($\\lambda={lam:.2f}$)')
    ax2.annotate('*', xy=(-w / 2, obs_hist[0] + 5), fontsize=14, ha='center',
                 color='#3498db', fontweight='bold')
    ax2.set_xlabel('Satellites per main star'); ax2.set_ylabel('Number of stars')
    ax2.set_xticks(x); ax2.set_xticklabels(k_range)
    all_c = [0] * N_zero + sats_per
    disp = np.var(all_c) / np.mean(all_c)
    ax2.set_title(f'(b) Poisson fit ($N = {N_true}$, disp. = {disp:.2f})', fontweight='bold')
    ax2.legend(fontsize=10); ax2.grid(True, alpha=0.15, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, 'p3_fig2.png'), dpi=300); plt.close()
    print("  Figure 2 saved.")


# ── Figure 3: Density vs n ──
def figure3():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    n_arr = np.array(stars); s_arr = np.array(sats_per)
    n_bins = np.linspace(50e9, 200e9, 31)
    n_cen = (n_bins[:-1] + n_bins[1:]) / 2
    means = [np.mean(s_arr[(n_arr >= n_bins[i]) & (n_arr < n_bins[i + 1])])
             if ((n_arr >= n_bins[i]) & (n_arr < n_bins[i + 1])).sum() > 0 else np.nan
             for i in range(len(n_bins) - 1)]
    ax1.scatter(n_arr / 1e9, s_arr, s=3, alpha=0.15, c='#3498db', rasterized=True)
    ax1.plot(n_cen / 1e9, means, 'o-', color='#e74c3c', markersize=4, linewidth=1.5,
             label='Binned mean')
    n_th = np.linspace(50e9, 200e9, 200)
    d_th = 46 * np.log10(n_th) + 1.67
    ax1.plot(n_th / 1e9, 5000 / (d_th * np.log(10)), '--', color='#27ae60', linewidth=1.5,
             label=r"Cram$\acute{\mathrm{e}}$r: $R/\ln P$")
    ax1.set_xlabel(r'$n$ (billions)'); ax1.set_ylabel('Satellites per star')
    ax1.set_title('(a) Satellite density vs $n$', fontweight='bold')
    ax1.legend(fontsize=9); ax1.set_ylim(0, 16); ax1.grid(True, alpha=0.2)

    n_bins2 = np.linspace(50e9, 200e9, 13)
    ratios, cen2 = [], []
    for i in range(len(n_bins2) - 1):
        mask = (n_arr >= n_bins2[i]) & (n_arr < n_bins2[i + 1])
        if mask.sum() < 10: continue
        d_m = np.mean(46 * np.log10(n_arr[mask]) + 1.67)
        ratios.append(np.mean(s_arr[mask]) / (5000 / (d_m * np.log(10))))
        cen2.append((n_bins2[i] + n_bins2[i + 1]) / 2)
    ax2.bar(range(len(ratios)), ratios, color='#9b59b6', alpha=0.7, edgecolor='#333')
    ax2.axhline(y=1.0, color='#e74c3c', linestyle='--', linewidth=1.5, label='Perfect agreement')
    ax2.set_xticks(range(len(ratios)))
    ax2.set_xticklabels([f'{c / 1e9:.0f}' for c in cen2], fontsize=8, rotation=45)
    ax2.set_xlabel('$n$ bin center (B)')
    ax2.set_ylabel("Observed / Cram\\'er")
    ax2.set_title("(b) Ratio by region", fontweight='bold')
    ax2.set_ylim(0.8, 1.3); ax2.legend(fontsize=10); ax2.grid(True, alpha=0.15, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, 'p3_fig3.png'), dpi=300); plt.close()
    print("  Figure 3 saved.")


# ── Figure 4: Close encounters + small gaps ──
def figure4():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    small = [(n, k) for n, k in satellites if k <= 100]
    ns, ks = zip(*small)
    ax1.scatter(np.array(ns) / 1e9, ks, c='gold', s=60, edgecolors='darkorange',
                linewidths=0.8, zorder=5, marker='*')
    for n, k in small:
        if k <= 20:
            ax1.annotate(f'k={k}', xy=(n / 1e9, k), fontsize=7,
                         xytext=(3, 3), textcoords='offset points')
    ax1.set_xlabel('Main star $n$ (billions)'); ax1.set_ylabel('Gap $k$')
    ax1.set_title('(a) Close encounters ($k \\leq 100$)', fontweight='bold')
    ax1.axhline(y=6, color='red', linestyle=':', alpha=0.4)
    ax1.text(55, 8, '$k=6$ (sexy primes)', fontsize=8, color='red')
    ax1.set_xlim(50, 205); ax1.grid(True, alpha=0.2)

    gap_counts = Counter(all_gaps)
    k_adm = [k for k in range(2, 62, 2) if k % 3 != 1]
    obs_c = [gap_counts.get(k, 0) for k in k_adm]
    ax2.bar(range(len(k_adm)), obs_c, color='#3498db', alpha=0.8, edgecolor='#333',
            linewidth=0.5)
    ax2.set_xticks(range(len(k_adm)))
    ax2.set_xticklabels([str(k) for k in k_adm], fontsize=7, rotation=45)
    ax2.set_xlabel('Gap $k$'); ax2.set_ylabel('Count (over 2107 stars)')
    ax2.set_title('(b) Fine-grained small-gap census ($k < 62$)', fontweight='bold')
    ax2.axhline(y=len(all_gaps) / 1667, color='grey', linestyle=':', alpha=0.5)
    ax2.grid(True, alpha=0.15, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, 'p3_fig4.png'), dpi=300); plt.close()
    print("  Figure 4 saved.")


if __name__ == '__main__':
    print("Generating figures for Paper IV (Satellite Primes)...")
    figure1()
    figure2()
    figure3()
    figure4()
    print("Done.")
