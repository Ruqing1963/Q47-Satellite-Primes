# Q47 Satellite Primes

**Part IV of the Titan Project**

## Satellite Primes: The Local Prime Landscape around Giant Primes from Q(n) = n⁴⁷ − (n−1)⁴⁷

A Cramér-Model Validation at 500-Digit Scales

**Author:** Ruqing Chen — GUT Geoservice Inc., Montréal, Canada

---

## Summary

This paper studies the local prime environment surrounding giant primes P = Q(n) = n⁴⁷ − (n−1)⁴⁷, probing the nearest primes within radius R = 5,000 of each P.

For **2,107** main-star primes of 494–521 digits, **9,012 satellite primes** P − k were discovered, matching the Cramér random model:

| Statistic | Cramér prediction | Observed | Ratio |
|:---|---:|---:|---:|
| Mean satellites/star | 4.35 | **4.28** | **0.98** |
| Poisson dispersion index | 1.000 | **1.07** | ~perfect |
| Gap uniformity χ² | — | p = 0.31 | not rejected |
| Nearest-satellite CDF | theory | observed | 1.01 |

---

## 🌟 Headline Result: Twin–Sexy Symmetry

A confirmation scan at R = 100 across all **2,992** main stars reveals:

| Gap k | Name | Count | E (conditional HL) | σ |
|:---|:---|---:|---:|---:|
| k = 2 | **Twin primes** | **7** | 7.2 | −0.1 |
| k = 6 | **Sexy primes** | **7** | 7.2 | −0.1 |
| k = 8 | Octet primes | 2 | 7.2 | −1.9 |

The equality **N_twin = N_sexy = 7** is a direct empirical confirmation of the identity **S_cond(k=2) ≡ S_cond(k=6)**, predicted by the Bayesian concentration principle: the fixed residue P ≡ 1 (mod 3) doubles the conditional twin-prime rate for k ≡ 2 (mod 6), exactly compensating the smaller unconditional singular series.

---

## 7 Twin Prime Pairs (500-digit scale)

| # | Main-star n | Approx. digits |
|---|---|---|
| 1 | 41,262,186,068 | ~498 |
| 2 | 63,150,957,871 | ~507 |
| 3 | 68,875,255,098 | ~509 |
| 4 | 123,037,305,946 | ~521 |
| 5 | 124,340,002,320 | ~521 |
| 6 | 126,720,185,653 | ~521 |
| 7 | 193,087,289,846 | ~530 |

## 7 Sexy Prime Pairs (500-digit scale)

| # | Main-star n | Approx. digits |
|---|---|---|
| 1 | 29,707,259,863 | ~492 |
| 2 | 103,957,400,503 | ~518 |
| 3 | 105,463,974,584 | ~518 |
| 4 | 122,726,858,404 | ~521 |
| 5 | 152,789,753,532 | ~524 |
| 6 | 154,849,622,427 | ~525 |
| 7 | 166,607,083,748 | ~526 |

---

## The Triple Coincidence: 7 = 7 = 7

Three independent phenomena in the Q(n) = n⁴⁷ − (n−1)⁴⁷ system all yield exactly 7:

| Phenomenon | Count | Source | Mechanism |
|---|---|---|---|
| Quintuplets | 7 | Part III | Bateman–Horn for k-tuples |
| Twin prime satellites | 7 | Part IV | Conditional Hardy–Littlewood |
| Sexy prime satellites | 7 | Part IV | Conditional Hardy–Littlewood |

The twin–sexy equality is theoretically predicted; the triple coincidence with quintuplets is numerical happenstance (P ≈ 1/308).

---

## The 3-Smooth Baseline Family

Gaps k = 2ᵃ × 3ᵇ (only prime factors 2 and 3) all share the **same** conditional singular series S_cond = 2.64 and the same expected count E ≈ 7.2. The 12 members within k ≤ 100:

| k | Factorization | Observed | σ |
|---|---|---|---|
| 2 | 2 | 7 ★ | −0.1 |
| 6 | 2×3 | 7 ★ | −0.1 |
| 8 | 2³ | 2 | −1.9 |
| 12 | 2²×3 | 6 | −0.4 |
| 18 | 2×3² | 7 ★ | −0.1 |
| 24 | 2³×3 | 3 | −1.6 |
| 32 | 2⁵ | 5 | −0.8 |
| 36 | 2²×3² | 6 | −0.4 |
| 48 | 2⁴×3 | 0 ⚠ | −2.7 |
| 54 | 2×3³ | 8 | +0.3 |
| 72 | 2³×3² | 8 | +0.3 |
| 96 | 2⁵×3 | 11 | +1.4 |

k = 8 having only 2 pairs (vs k = 2 having 7) is pure Poisson fluctuation — both have identical E = 7.2.

---

## Key Theoretical Results

### Forbidden Residue Lattice
All gaps satisfy k ≡ 0 or 2 (mod 6), because Q(n) ≡ 1 (mod 6) for all n. This eliminates 1/3 of even gaps, concentrating satellite density to ~3/ln(P) per admissible slot.

### Bayesian Concentration Principle
- **k ≡ 2 (mod 6)**: Bayesian factor B = 2 (all such pairs require P ≡ 1 mod 3)
- **k ≡ 0 (mod 6)**: Bayesian factor B = 1 (no concentration)
- The doubling for k ≡ 2 exactly compensates its smaller unconditional HL factor
- Result: mod-6 satellite classes are equal (4,468 vs 4,544; ratio 1.02)

---

## Repository Structure

```
Q47-Satellite-Primes/
├── paper/
│   ├── Q47_Satellite_Primes.tex        # LaTeX source (14 pages)
│   ├── Q47_Satellite_Primes.pdf        # Compiled paper
│   └── figures/
│       ├── p3_fig1.{pdf,png}           # Gap distribution + mod-30 structure
│       ├── p3_fig2.{pdf,png}           # Nearest-satellite CDF + Poisson fit
│       ├── p3_fig3.{pdf,png}           # Density vs main-star size
│       └── p3_fig4.{pdf,png}           # Close encounters detail
├── scripts/
│   ├── titan_radar_ultimate_5000.py    # R=5000 deep scan (satellite discovery)
│   ├── titan_radar_ultimate_100.py     # R=100 confirmation scan (twin/sexy census)
│   ├── analyze_satellites.py           # Statistical analysis pipeline
│   └── generate_figures.py             # Figure generation
├── data/
│   ├── satellites_9012.csv             # All 9,012 satellites (R=5000, 2,079 stars)
│   ├── star_summary_2079.csv           # Per-star satellite counts
│   ├── close_encounters_r100.csv       # All gaps k≤100 (R=100, 2,992 stars)
│   ├── twin_primes_7.csv              # Complete twin prime catalog
│   ├── sexy_primes_7.csv              # Complete sexy prime catalog
│   ├── smooth_baseline_family.csv      # 3-smooth baseline analysis
│   ├── conditional_hardy_littlewood.csv # Conditional HL table (k≤30)
│   ├── poisson_fit.csv                 # Poisson fit data
│   ├── density_by_range.csv            # Density by n-range
│   └── mod30_distribution.csv          # Mod-30 gap distribution
├── README.md
├── LICENSE
└── .gitignore
```

---

## Data Sources

| Dataset | Stars | Radius | Satellites | Purpose |
|---|---|---|---|---|
| R = 5000 (partial) | 2,079 (+28 zero-sat) | [2, 5000] | 9,012 | Cramér model validation |
| R = 100 (complete) | 2,992 | [2, 100] | 235 | Twin/sexy prime census |

The R = 5000 scan covered n ∈ [5.29×10¹⁰, 2.00×10¹¹] (2,107 stars inferred).
The R = 100 scan covered the full quadruplet catalog: 748 × 4 = 2,992 stars, n ∈ [2.19×10⁸, 2.00×10¹¹].

---

## Reproducing the Results

Each twin/sexy prime can be verified independently:
```python
from sympy import isprime

def Q(n):
    return n**47 - (n-1)**47

# Example: first twin prime pair
n = 41262186068
P = Q(n)
print(f"P has {len(str(P))} digits")
print(f"P is prime: {isprime(P)}")
print(f"P-2 is prime: {isprime(P-2)}")  # Twin!
```

---

## Titan Project Series

| Part | Title | Status | DOI |
|---|---|---|---|
| I | Statistical Morphology | Published | [10.5281/zenodo.18701355](https://zenodo.org/records/18701355) |
| II | Quadruplet Census | Published | [10.5281/zenodo.18728540](https://zenodo.org/records/18728540) |
| III | Quintuplet Boundary | Published | [10.5281/zenodo.18728917](https://zenodo.org/records/18728917) |
| **IV** | **Satellite Primes** | **Ready for upload** | *Pending* |
| V | Deep Sieve Structure | Planned | — |

---

## License

MIT License. See [LICENSE](LICENSE).
