# BenfordMA Design Spec — Benford's Law Fabrication Screening for Meta-Analysis

**Date:** 2026-03-25
**Target:** `C:\BenfordMA\benford-ma.html` (new single-file HTML app)
**Build step:** `C:\BenfordMA\build_digits.py` (new Python script)
**Data sources:** FragilityAtlas + PredictionGap study-level data (403 Cochrane reviews)

## Background

Benford's Law states that the first significant digit of naturally occurring numbers follows `P(d) = log10(1 + 1/d)` for d=1..9. Data that spans multiple orders of magnitude (effect sizes, SEs, sample sizes) should conform. Fabricated data deviates because humans are poor at generating "natural-looking" digit distributions.

This has been applied extensively in forensic accounting and election fraud detection, but never systematically to meta-analysis at scale. We apply it to 403 Cochrane systematic reviews to screen for potential data fabrication.

## 1. Data Pipeline

### `build_digits.py`

Extract first and second significant digits from 4 numeric fields across all studies in the 403 reviews.

**Fields:**
1. Effect sizes (theta from FragilityAtlas: `median_theta` or study-level effects)
2. Standard errors (derived from CI width or directly available)
3. Sample sizes / study k (number of participants or studies)
4. P-values (from PredictionGap: `p_value`, or computed from effect/SE as `2*(1 - normalCDF(|effect/SE|))`)

**Source files:**
- `C:\FragilityAtlas\data\output\fragility_atlas_results.csv` (403 reviews, review-level aggregates)
- `C:\PredictionGap\data\output\prediction_gap_results.csv` (403 reviews, includes theta, CI, p-value, tau2, I2)

**Note:** These files contain review-level aggregates, not individual study data. Each review contributes one set of values (pooled effect, pooled SE, k, p-value). With 403 reviews x 4 fields = up to 1,612 digit observations for corpus-level analysis. For per-review analysis, we analyze the distribution of digits across the 4 fields within that review (low power, supplemented by Fisher's method across reviews).

**If study-level data is available** (check for `C:\FragilityAtlas\data\output\study_level_*.csv` or similar), use that instead for much higher power.

**Digit extraction:**
```python
def first_digit(x):
    """Return first significant digit (1-9) of abs(x), or None."""
    if x is None or x == 0:
        return None
    x = abs(x)
    while x < 1:
        x *= 10
    while x >= 10:
        x /= 10
    return int(x)

def second_digit(x):
    """Return second significant digit (0-9) of abs(x), or None."""
    if x is None or x == 0:
        return None
    x = abs(x)
    while x < 1:
        x *= 10
    while x >= 10:
        x /= 10
    return int((x * 10) % 10)
```

**Output:** `data/corpus_digits.json`

```json
{
  "reviews": [
    {
      "review_id": "CD000028",
      "analysis_name": "...",
      "k": 21,
      "digits": [
        {"field": "effect", "value": -0.2945, "d1": 2, "d2": 9},
        {"field": "se", "value": 0.0741, "d1": 7, "d2": 4},
        {"field": "n", "value": 21, "d1": 2, "d2": 1},
        {"field": "pvalue", "value": 0.000003, "d1": 3, "d2": 0}
      ]
    }
  ],
  "corpus_summary": {
    "n_reviews": 403,
    "n_digits": 1612,
    "fields": ["effect", "se", "n", "pvalue"]
  }
}
```

## 2. Statistical Methods

### Benford's expected distribution

First digit: `P(d) = log10(1 + 1/d)` for d = 1..9
```
d:  1     2     3     4     5     6     7     8     9
P:  30.1% 17.6% 12.5% 9.7%  7.9%  6.7%  5.8%  5.1%  4.6%
```

Second digit: `P(d) = sum_{k=1..9} log10(1 + 1/(10k+d))` for d = 0..9
```
d:  0     1     2     3     4     5     6     7     8     9
P:  12.0% 11.4% 10.9% 10.4% 10.0% 9.7%  9.3%  9.0%  8.8%  8.5%
```

### Corpus-level tests

1. **Chi-squared goodness-of-fit**: `chi2 = sum((O_d - E_d)^2 / E_d)` with df=8 (first digit) or df=9 (second digit). Report chi2, p-value.

2. **Mean Absolute Deviation (MAD)**: `MAD = (1/K) * sum(|O_d/N - P_d|)` where K=9 (first) or K=10 (second).
   - Nigrini classification: Conforming (<0.006), Acceptable (0.006-0.012), Marginally acceptable (0.012-0.015), Nonconforming (>0.015)

3. **Per-field breakdown**: Separate chi2 and MAD for each of the 4 fields. Effect sizes and SEs should conform well (span orders of magnitude). P-values may deviate legitimately (bounded 0-1, cluster near thresholds). Sample sizes may deviate if studies cluster in narrow ranges.

4. **Mantissa Arc Test**: For each value x, compute `m = log10(|x|) mod 1`. Under Benford, mantissas are uniformly distributed on [0,1). Test uniformity via Kolmogorov-Smirnov. Visualize on polar plot (uniform = circle, fabrication = spikes).

### Per-review screening

Each review has ~4 digit observations (one per field), which is too few for chi-squared. Instead:

1. **Exact binomial test per digit**: For each first digit d, test whether this review's proportion of d differs from Benford P(d). With only 4 observations per review, this has very low power individually.

2. **Fisher's method across fields**: Combine p-values from binomial tests using Fisher's method: `X = -2 * sum(ln(p_i))`, which follows chi2(2k) under the null.

3. **FDR correction**: Apply Benjamini-Hochberg across 403 reviews. Flag reviews with FDR-adjusted p < 0.05.

4. **Suspicion score (0-100)**:
   - MAD contribution (0-40): `min(40, MAD / 0.025 * 40)` — maps MAD to [0, 0.025] range
   - Fisher combined p contribution (0-40): `min(40, -log10(fisher_p) / 4 * 40)` — maps p from 1 to 1e-4
   - Caliper ratio contribution (0-20): if p-value in [0.04, 0.05), add 20; if in [0.03, 0.04), add 10; else 0

5. **Classification:**
   - Conforming: score < 25
   - Minor Deviation: score 25-50
   - Suspicious: score 50-75
   - Highly Suspicious: score > 75

**Important caveat:** Low per-review power means most reviews will be "Conforming" even if mildly fabricated. The tool is a screening flag, not proof of fabrication. This must be stated prominently.

### P-value caliper test (complementary)

Not strictly Benford, but a well-established fabrication signal:
- Count p-values in [0.04, 0.05) vs [0.05, 0.06)
- Under the null (no fabrication), these bins should have similar counts
- Excess in [0.04, 0.05) suggests "nudging" results to significance
- Test: exact binomial comparing the two bins

## 3. Dashboard Layout

Single scrollable page, dark mode toggle. Two modes: Corpus (embedded data) and Custom (user paste).

### Header

- Title: "BenfordMA — Digit Forensics for Meta-Analysis"
- Subtitle: "Benford's Law screening of N Cochrane reviews"
- Mode toggle: [Corpus | Custom]

### Corpus Mode

**Section 1 — Corpus Summary**: 4 stat boxes:
1. Total reviews (403)
2. Total digit observations (N)
3. Overall first-digit MAD (with Nigrini classification badge)
4. Overall chi-squared p-value

**Section 2 — First Digit Distribution**:
- Main chart: 9 bars (observed) with Benford expected line overlay (SVG)
- Chi-squared and MAD below the chart
- 2x2 grid of per-field charts (Effect, SE, N, P-value), each smaller

**Section 3 — Second Digit Distribution**: Same layout as Section 2 for second digits (10 bars).

**Section 4 — Mantissa Arc Plot**: Polar/circular histogram of mantissa values. SVG circle with 36 bins (10-degree sectors). Uniform = even circle, fabrication = spikes.

**Section 5 — P-value Caliper**: Histogram of p-values in bins of width 0.01 from 0 to 0.10. Highlight [0.04, 0.05) and [0.05, 0.06) bins. Caliper test result.

**Section 6 — Per-Review Screening Table**: Sortable table:
| Review ID | Analysis | k | MAD | Fisher p | FDR p | Score | Classification |
Click to expand: that review's first-digit bar chart vs Benford expected.

### Custom Mode

- CSV paste area: "Study, Effect, SE, N, P-value"
- Import button
- Same Sections 1-5 render on the pasted data
- No Section 6 (single analysis, not per-review)

### Footer

- "BenfordMA v1.0 — Browser-based, no data leaves your device."
- Export: [CSV] filtered table, [PNG] first-digit chart
- Caveat text: "Benford deviation is a screening signal, not proof of fabrication. Low digit counts reduce statistical power."

## 4. Visual Design

Same CSS custom properties pattern as other portfolio apps:
- Light/dark theming with `data-theme` attribute
- `eqd_` → `bma_` localStorage prefix
- Card-based layout, responsive grid
- Print styles hiding interactive elements

## 5. Integration Map

| File | Purpose |
|------|---------|
| `C:\BenfordMA\build_digits.py` | Extract digits from FragilityAtlas + PredictionGap data |
| `C:\BenfordMA\data\corpus_digits.json` | Build artifact: digit data for 403 reviews |
| `C:\BenfordMA\benford-ma.html` | Single-file HTML dashboard (~2,000 lines) |

## 6. Out of Scope

- Bayesian digit analysis (Pericchi-Torres)
- Network analysis of fabrication patterns across reviews
- Integration with COPE or Retraction Watch databases
- Automated remediation or study removal recommendations
- Second-order tests (summation, distortion factor)

## 7. Validation

- Verify Benford expected frequencies sum to 1.0
- Simulate 10,000 log-normal random numbers: confirm first-digit MAD < 0.006
- Simulate uniform random digits 1-9: confirm MAD > 0.015 (nonconforming)
- Chi-squared df correct: 8 for first digit, 9 for second digit
- Fisher's method df correct: 2k where k = number of tests combined
- P-value caliper: equal-sized bins, test is two-sided binomial
- Per-review classification thresholds produce reasonable distribution (not all Conforming, not all Suspicious)
