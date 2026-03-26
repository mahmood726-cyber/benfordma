# BenfordMA Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a single-file HTML dashboard that applies Benford's Law digit analysis to 403 Cochrane meta-analyses, screening for potential data fabrication at corpus and per-review levels.

**Architecture:** Two-phase build: (1) Python script extracts first/second digits from review-level aggregates (403 reviews x 6 numeric fields = ~2,400 digit observations) plus 394K specification-level values as sensitivity analysis. (2) Single-file HTML dashboard with Benford distribution charts, MAD tests, mantissa arc plot, p-value caliper, and per-review screening table.

**Tech Stack:** Python 3.x (csv, json, math), vanilla HTML/CSS/JS, SVG for charts.

**Spec:** `C:\BenfordMA\docs\superpowers\specs\2026-03-25-benford-ma-design.md`

---

## File Structure

| File | Responsibility |
|------|---------------|
| `C:\BenfordMA\build_digits.py` | Extract digits from FragilityAtlas + PredictionGap data |
| `C:\BenfordMA\data\corpus_digits.json` | Build artifact: digit data for 403 reviews |
| `C:\BenfordMA\benford-ma.html` | Single-file HTML dashboard (~2,000 lines) |

## Data Sources

| Source | Path | Rows | Type | Fields used |
|--------|------|------|------|-------------|
| FA results | `C:\FragilityAtlas\data\output\fragility_atlas_results.csv` | 403 | Review-level | median_theta, iqr_theta, k, frac_significant |
| PG results | `C:\PredictionGap\data\output\prediction_gap_results.csv` | 403 | Review-level | theta, ci_lo, ci_hi, tau2, I2, p_value |
| FA specs | `C:\FragilityAtlas\data\output\fragility_atlas_specifications.csv` | 394,570 | Spec-level | theta, se_theta, p_value (sensitivity analysis) |
| FA validation | `C:\FragilityAtlas\data\output\r_validation_inputs.json` | 10 analyses | Per-study | yi[], sei[] arrays |

---

### Task 1: Build `build_digits.py` — extract digit data

**Files:**
- Create: `C:\BenfordMA\build_digits.py`
- Output: `C:\BenfordMA\data\corpus_digits.json`

**Context:** Extract first and second significant digits from numeric fields across 403 reviews. Primary source is the review-level CSVs (theta, SE, p-value, tau2, I2, k). Also extract from the 394K specification file for a sensitivity analysis corpus. Per-study data from r_validation_inputs.json provides a ground-truth validation subset.

- [ ] **Step 1: Write build_digits.py**

```python
"""Extract digit data for Benford's Law analysis of Cochrane meta-analyses.

Sources:
  - FragilityAtlas results (403 reviews, 6 numeric fields)
  - PredictionGap results (403 reviews, theta + CI + tau2 + I2 + p)
  - FragilityAtlas specifications (394K rows, sensitivity analysis)
  - r_validation_inputs.json (10 analyses, per-study yi/sei)

Output: data/corpus_digits.json
"""

import csv
import json
import math
from pathlib import Path


def first_digit(x):
    """Return first significant digit (1-9) of abs(x), or None."""
    if x is None or not math.isfinite(x) or x == 0:
        return None
    x = abs(x)
    while x < 1:
        x *= 10
    while x >= 10:
        x /= 10
    return int(x)


def second_digit(x):
    """Return second significant digit (0-9) of abs(x), or None."""
    if x is None or not math.isfinite(x) or x == 0:
        return None
    x = abs(x)
    while x < 1:
        x *= 10
    while x >= 10:
        x /= 10
    return int((x * 10) % 10)


def mantissa(x):
    """Return mantissa = log10(|x|) mod 1, or None."""
    if x is None or not math.isfinite(x) or x == 0:
        return None
    return math.log10(abs(x)) % 1


def safe_float(val):
    if val is None or val == '':
        return None
    try:
        v = float(val)
        return v if math.isfinite(v) else None
    except (ValueError, TypeError):
        return None


def load_csv(path):
    rows = []
    p = Path(path)
    if not p.exists():
        print(f"  WARNING: {path} not found")
        return rows
    with open(p, encoding='utf-8') as f:
        for row in csv.DictReader(f):
            rows.append(row)
    print(f"  Loaded {len(rows)} rows from {p.name}")
    return rows


def extract_review_digits(fa_rows, pg_rows):
    """Extract digits from review-level data. Returns list of review records."""
    # Index PG by review_id
    pg_map = {}
    for r in pg_rows:
        pg_map[r['review_id']] = r

    reviews = []
    for fa in fa_rows:
        rid = fa['review_id']
        pg = pg_map.get(rid, {})

        digits = []
        # Fields to extract from both sources
        fields = {
            'effect': safe_float(pg.get('theta')) or safe_float(fa.get('median_theta')),
            'se': safe_float(fa.get('iqr_theta')),  # proxy for SE spread
            'k': safe_float(fa.get('k')),
            'pvalue': safe_float(pg.get('p_value')),
            'tau2': safe_float(pg.get('tau2')),
            'I2': safe_float(pg.get('I2')),
        }

        for field_name, val in fields.items():
            if val is not None and val != 0:
                d1 = first_digit(val)
                d2 = second_digit(val)
                m = mantissa(val)
                if d1 is not None:
                    digits.append({
                        'field': field_name,
                        'value': round(val, 8),
                        'd1': d1,
                        'd2': d2,
                        'mantissa': round(m, 6) if m is not None else None,
                    })

        reviews.append({
            'review_id': rid,
            'analysis_name': fa.get('analysis_name', ''),
            'k': int(float(fa.get('k', 0))),
            'digits': digits,
        })

    return reviews


def extract_spec_digits(spec_path):
    """Extract digits from 394K specification-level data (sensitivity analysis)."""
    specs = load_csv(spec_path)
    d1_counts = [0] * 10  # d1_counts[1..9]
    d2_counts = [0] * 10  # d2_counts[0..9]
    n_total = 0

    for row in specs:
        for col in ['theta', 'se_theta', 'p_value']:
            val = safe_float(row.get(col))
            if val is not None and val != 0:
                d1 = first_digit(val)
                d2 = second_digit(val)
                if d1 is not None:
                    d1_counts[d1] += 1
                    n_total += 1
                if d2 is not None:
                    d2_counts[d2] += 1

    return {
        'n_values': n_total,
        'd1_counts': d1_counts[1:],  # [count_for_1, ..., count_for_9]
        'd2_counts': d2_counts,       # [count_for_0, ..., count_for_9]
    }


def extract_study_digits(json_path):
    """Extract per-study digits from r_validation_inputs.json (ground truth)."""
    p = Path(json_path)
    if not p.exists():
        print(f"  WARNING: {json_path} not found")
        return None

    with open(p, encoding='utf-8') as f:
        data = json.load(f)

    study_digits = []
    for analysis in data:
        yi = analysis.get('yi', [])
        sei = analysis.get('sei', [])
        for i, (y, s) in enumerate(zip(yi, sei)):
            for field, val in [('effect', y), ('se', s)]:
                if val is not None and val != 0:
                    d1 = first_digit(val)
                    d2 = second_digit(val)
                    m = mantissa(val)
                    if d1 is not None:
                        study_digits.append({
                            'review_id': analysis.get('review_id', ''),
                            'study_idx': i,
                            'field': field,
                            'value': round(val, 8),
                            'd1': d1,
                            'd2': d2,
                            'mantissa': round(m, 6) if m is not None else None,
                        })

    print(f"  Extracted {len(study_digits)} per-study digit observations")
    return study_digits


def main():
    print("Building BenfordMA digit corpus...\n")

    # Review-level digits (primary)
    print("1. Review-level data:")
    fa = load_csv(r'C:\FragilityAtlas\data\output\fragility_atlas_results.csv')
    pg = load_csv(r'C:\PredictionGap\data\output\prediction_gap_results.csv')
    reviews = extract_review_digits(fa, pg)

    total_digits = sum(len(r['digits']) for r in reviews)
    print(f"  -> {len(reviews)} reviews, {total_digits} digit observations\n")

    # Specification-level sensitivity (secondary)
    print("2. Specification-level sensitivity analysis:")
    spec_summary = extract_spec_digits(
        r'C:\FragilityAtlas\data\output\fragility_atlas_specifications.csv'
    )
    print(f"  -> {spec_summary['n_values']} specification-level digit observations\n")

    # Per-study ground truth (tertiary)
    print("3. Per-study ground truth:")
    study_digits = extract_study_digits(
        r'C:\FragilityAtlas\data\output\r_validation_inputs.json'
    )

    # Build output
    output = {
        'reviews': reviews,
        'spec_sensitivity': spec_summary,
        'study_ground_truth': study_digits,
        'corpus_summary': {
            'n_reviews': len(reviews),
            'n_review_digits': total_digits,
            'n_spec_digits': spec_summary['n_values'],
            'n_study_digits': len(study_digits) if study_digits else 0,
            'fields': ['effect', 'se', 'k', 'pvalue', 'tau2', 'I2'],
        }
    }

    out_path = Path(r'C:\BenfordMA\data\corpus_digits.json')
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(output, f)

    print(f"\n{'='*50}")
    print("BENFORD DIGIT CORPUS BUILT")
    print(f"{'='*50}")
    print(f"  Reviews: {len(reviews)}")
    print(f"  Review-level digits: {total_digits}")
    print(f"  Spec-level digits: {spec_summary['n_values']}")
    print(f"  Per-study digits: {len(study_digits) if study_digits else 0}")
    print(f"  Output: {out_path}")
    print(f"  Size: {out_path.stat().st_size / 1024:.0f} KB")


if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Run build script**

```bash
cd /c/BenfordMA && python build_digits.py
```

Expected: 403 reviews, ~2,000+ review-level digits, ~1M spec-level digits, ~640 per-study digits.

- [ ] **Step 3: Verify JSON structure**

```bash
python -c "import json; d=json.load(open('data/corpus_digits.json')); print(f'Reviews: {d[\"corpus_summary\"][\"n_reviews\"]}'); print(f'Review digits: {d[\"corpus_summary\"][\"n_review_digits\"]}'); print(f'Spec digits: {d[\"corpus_summary\"][\"n_spec_digits\"]}'); print(json.dumps(d['reviews'][0], indent=2)[:500])"
```

- [ ] **Step 4: Commit**

```bash
cd /c/BenfordMA && git add build_digits.py data/corpus_digits.json && git commit -m "feat: build digit extraction pipeline — 403 reviews + 394K specs + per-study"
```

---

### Task 2: Dashboard scaffold — HTML, CSS, dark mode, embedded data

**Files:**
- Create: `C:\BenfordMA\benford-ma.html`

**Context:** Create the full HTML shell with all sections, CSS theming, mode toggle (Corpus/Custom), embedded JSON data. Follow the same single-file pattern as RMST Meta and EvidenceQuality dashboards.

- [ ] **Step 1: Create benford-ma.html**

CSS requirements:
- Custom properties for light/dark theming (`bma_` localStorage prefix)
- `.stat-grid` (4 stat boxes), `.card`, `.btn`
- `.benford-bar` class for chart bars (observed vs expected)
- `.grade-conforming` (green), `.grade-acceptable` (blue), `.grade-marginal` (amber), `.grade-nonconforming` (red)
- `.sortable-table`, `.detail-panel` for accordion
- Print styles, responsive at 768px

HTML structure:
- Header: "BenfordMA — Digit Forensics for Meta-Analysis"
- Mode toggle: [Corpus | Custom] radio buttons
- Corpus mode: Sections 1-6 (summary, first digit, second digit, mantissa, caliper, table)
- Custom mode: CSV paste area + import button, same sections 1-5
- Footer: "BenfordMA v1.0 — Browser-based, no data leaves your device."
- Caveat text: "Benford deviation is a screening signal, not proof of fabrication."
- Export buttons: CSV + PNG

JS scaffold:
- `var DATA = <embedded JSON>;`
- `var filtered = [];`
- `escapeHtml()`, `applyMode()`, dark mode toggle
- Stubs for all render functions

- [ ] **Step 2: Embed actual JSON data from `data/corpus_digits.json`**

Read the JSON file and paste as `var DATA = <json>;` in the script block. Verify no `</script>` in the data.

- [ ] **Step 3: Verify — page loads, mode toggle works, dark mode toggles, no console errors**

- [ ] **Step 4: Commit**

```bash
cd /c/BenfordMA && git add benford-ma.html && git commit -m "feat: dashboard scaffold — HTML, CSS, dark mode, mode toggle, embedded data"
```

---

### Task 3: Benford statistics engine + summary cards

**Files:**
- Modify: `C:\BenfordMA\benford-ma.html` — JS section

**Context:** Implement all the statistical functions and the corpus summary. These are pure functions that take digit arrays and return test results.

- [ ] **Step 1: Implement Benford statistical functions**

```javascript
/* Benford expected frequencies */
function benfordFirst() {
  // P(d) = log10(1 + 1/d) for d=1..9
  var p = [];
  for (var d = 1; d <= 9; d++) p.push(Math.log10(1 + 1 / d));
  return p;
}

function benfordSecond() {
  // P(d) = sum_{k=1..9} log10(1 + 1/(10k+d)) for d=0..9
  var p = [];
  for (var d = 0; d <= 9; d++) {
    var s = 0;
    for (var k = 1; k <= 9; k++) s += Math.log10(1 + 1 / (10 * k + d));
    p.push(s);
  }
  return p;
}

/* Chi-squared goodness of fit */
function chiSquaredTest(observed, expected, n) {
  var chi2 = 0;
  for (var i = 0; i < observed.length; i++) {
    var e = expected[i] * n;
    if (e > 0) chi2 += (observed[i] - e) * (observed[i] - e) / e;
  }
  var df = observed.length - 1;
  var p = 1 - chi2CDF(chi2, df);
  return { chi2: chi2, df: df, p: p };
}

/* Mean Absolute Deviation */
function computeMAD(observed, expected, n) {
  var mad = 0;
  for (var i = 0; i < observed.length; i++) {
    mad += Math.abs(observed[i] / n - expected[i]);
  }
  return mad / observed.length;
}

/* Nigrini classification */
function nigriniClass(mad) {
  if (mad < 0.006) return 'Conforming';
  if (mad < 0.012) return 'Acceptable';
  if (mad < 0.015) return 'Marginally Acceptable';
  return 'Nonconforming';
}

/* chi2CDF — same as RMST Meta */
function chi2CDF(x, df) {
  if (x <= 0) return 0;
  // Use regularized gamma function approximation
  var k = df / 2, z = x / 2;
  // Series expansion for lower incomplete gamma
  var sum = 0, term = 1 / k;
  sum = term;
  for (var i = 1; i < 200; i++) {
    term *= z / (k + i);
    sum += term;
    if (Math.abs(term) < 1e-12) break;
  }
  return sum * Math.exp(-z + k * Math.log(z) - logGamma(k));
}

function logGamma(x) {
  // Stirling approximation
  var c = [76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  var y = x, tmp = x + 5.5;
  tmp -= (x + 0.5) * Math.log(tmp);
  var ser = 1.000000000190015;
  for (var j = 0; j < 6; j++) ser += c[j] / ++y;
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}

/* Kolmogorov-Smirnov test for mantissa uniformity */
function ksMantissa(mantissas) {
  var sorted = mantissas.slice().sort(function(a, b) { return a - b; });
  var n = sorted.length;
  var dMax = 0;
  for (var i = 0; i < n; i++) {
    var d1 = Math.abs((i + 1) / n - sorted[i]);
    var d2 = Math.abs(sorted[i] - i / n);
    dMax = Math.max(dMax, d1, d2);
  }
  // Approximate p-value
  var sqrtN = Math.sqrt(n);
  var lambda = (sqrtN + 0.12 + 0.11 / sqrtN) * dMax;
  var p = 2 * Math.exp(-2 * lambda * lambda);
  return { D: dMax, p: Math.min(1, Math.max(0, p)) };
}

/* P-value caliper test */
function caliperTest(pvalues) {
  var bin04 = 0, bin05 = 0;
  pvalues.forEach(function(p) {
    if (p >= 0.04 && p < 0.05) bin04++;
    if (p >= 0.05 && p < 0.06) bin05++;
  });
  // Exact binomial: under null, each falls in either bin with equal probability
  var total = bin04 + bin05;
  if (total === 0) return { bin04: 0, bin05: 0, ratio: null, p: 1 };
  // Binomial p-value: P(X >= bin04) under p=0.5
  var p = 0;
  for (var x = Math.max(bin04, bin05); x <= total; x++) {
    p += binomPMF(total, x, 0.5);
  }
  p = Math.min(1, 2 * p); // two-sided
  return {
    bin04: bin04,
    bin05: bin05,
    ratio: bin05 > 0 ? bin04 / bin05 : null,
    total: total,
    p: p,
  };
}

function binomPMF(n, k, p) {
  return Math.exp(logChoose(n, k) + k * Math.log(p) + (n - k) * Math.log(1 - p));
}

function logChoose(n, k) {
  return logGamma(n + 1) - logGamma(k + 1) - logGamma(n - k + 1);
}

/* Fisher's method for combining p-values */
function fisherCombine(pvalues) {
  var valid = pvalues.filter(function(p) { return p > 0 && p <= 1; });
  if (valid.length === 0) return { chi2: 0, df: 0, p: 1 };
  var chi2 = 0;
  valid.forEach(function(p) { chi2 += -2 * Math.log(p); });
  var df = 2 * valid.length;
  var p = 1 - chi2CDF(chi2, df);
  return { chi2: chi2, df: df, p: p };
}
```

- [ ] **Step 2: Implement `computeCorpusStats()` and `renderSummary()`**

```javascript
function computeCorpusStats() {
  // Collect all digits from reviews
  var allD1 = [], allD2 = [], allMantissas = [], allPvalues = [];
  var fieldD1 = { effect: [], se: [], k: [], pvalue: [], tau2: [], I2: [] };

  DATA.reviews.forEach(function(r) {
    r.digits.forEach(function(d) {
      if (d.d1 !== null) allD1.push(d.d1);
      if (d.d2 !== null) allD2.push(d.d2);
      if (d.mantissa !== null) allMantissas.push(d.mantissa);
      if (d.field === 'pvalue' && d.value !== null) allPvalues.push(d.value);
      if (fieldD1[d.field] && d.d1 !== null) fieldD1[d.field].push(d.d1);
    });
  });

  // Count first digits
  var d1Counts = new Array(9).fill(0);
  allD1.forEach(function(d) { d1Counts[d - 1]++; });
  var d2Counts = new Array(10).fill(0);
  allD2.forEach(function(d) { d2Counts[d]++; });

  var bfFirst = benfordFirst();
  var bfSecond = benfordSecond();
  var n1 = allD1.length;
  var n2 = allD2.length;

  return {
    n: n1,
    d1Counts: d1Counts,
    d2Counts: d2Counts,
    chi2First: chiSquaredTest(d1Counts, bfFirst, n1),
    chi2Second: chiSquaredTest(d2Counts, bfSecond, n2),
    madFirst: computeMAD(d1Counts, bfFirst, n1),
    madSecond: computeMAD(d2Counts, bfSecond, n2),
    mantissa: ksMantissa(allMantissas),
    caliper: caliperTest(allPvalues),
    fieldD1: fieldD1,
    allMantissas: allMantissas,
  };
}

function renderSummary(stats) {
  document.getElementById('statTotal').textContent = DATA.corpus_summary.n_reviews;
  document.getElementById('statDigits').textContent = stats.n;
  document.getElementById('statMAD').textContent = stats.madFirst.toFixed(4);
  document.getElementById('statClass').textContent = nigriniClass(stats.madFirst);
  // Color the classification badge
  var cls = nigriniClass(stats.madFirst);
  var badge = document.getElementById('statClass');
  badge.className = 'stat-value ' +
    (cls === 'Conforming' ? 'grade-conforming' :
     cls === 'Acceptable' ? 'grade-acceptable' :
     cls === 'Marginally Acceptable' ? 'grade-marginal' : 'grade-nonconforming');
}
```

Call `computeCorpusStats()` on load, store result, pass to all render functions.

- [ ] **Step 3: Verify — summary cards show correct totals and MAD**

- [ ] **Step 4: Commit**

```bash
cd /c/BenfordMA && git add benford-ma.html && git commit -m "feat: Benford statistics engine + summary cards"
```

---

### Task 4: First digit + second digit distribution charts

**Files:**
- Modify: `C:\BenfordMA\benford-ma.html` — `renderFirstDigit()`, `renderSecondDigit()`

**Context:** SVG bar charts showing observed vs Benford expected. Main chart for all digits, then 2x2 grid for per-field breakdown. Chi-squared and MAD displayed below each chart.

- [ ] **Step 1: Implement `renderFirstDigit(stats)`**

SVG bar chart: 9 bars (digits 1-9). Observed as filled bars, Benford expected as line/dots overlay. Width ~600px. Colors: bars use accent color, expected line in red/orange.

Below the main chart: chi-squared statistic, p-value, MAD, Nigrini class.

Per-field 2x2 grid: smaller charts (250px wide) for effect, SE, k, pvalue fields. Each with its own chi2/MAD.

- [ ] **Step 2: Implement `renderSecondDigit(stats)`**

Same layout but 10 bars (digits 0-9) with second-digit Benford expected.

- [ ] **Step 3: Verify — charts visible, Benford expected line shows characteristic 30.1% → 4.6% decline**

- [ ] **Step 4: Commit**

```bash
cd /c/BenfordMA && git add benford-ma.html && git commit -m "feat: first and second digit distribution charts with per-field breakdown"
```

---

### Task 5: Mantissa arc plot + p-value caliper

**Files:**
- Modify: `C:\BenfordMA\benford-ma.html` — `renderMantissa()`, `renderCaliper()`

**Context:** Mantissa arc is a polar histogram of log10(|x|) mod 1. P-value caliper is a histogram of p-values near 0.05 with the caliper test result.

- [ ] **Step 1: Implement `renderMantissa(stats)`**

SVG polar plot: circle divided into 36 sectors (10 degrees each). Radius of each sector proportional to count of mantissas in that bin. Uniform distribution = perfect circle. Draw a reference circle for the expected uniform count.

KS test result displayed below: D statistic, p-value, interpretation.

- [ ] **Step 2: Implement `renderCaliper(stats)`**

SVG histogram of p-values from 0 to 0.10 in bins of 0.01. Highlight [0.04, 0.05) in amber and [0.05, 0.06) in blue. Show caliper test: ratio, binomial p-value, interpretation.

- [ ] **Step 3: Verify — mantissa plot is roughly circular, caliper bins colored correctly**

- [ ] **Step 4: Commit**

```bash
cd /c/BenfordMA && git add benford-ma.html && git commit -m "feat: mantissa arc plot + p-value caliper histogram"
```

---

### Task 6: Per-review screening table + accordion

**Files:**
- Modify: `C:\BenfordMA\benford-ma.html` — `renderReviewTable()`, `toggleReviewDetail()`

**Context:** Sortable table of 403 reviews with suspicion scores. Click to expand showing that review's digit distribution. BH-FDR correction across 403 reviews.

- [ ] **Step 1: Compute per-review scores**

```javascript
function computePerReviewScores() {
  var bfFirst = benfordFirst();
  var results = DATA.reviews.map(function(r) {
    var d1s = r.digits.filter(function(d) { return d.d1 !== null; }).map(function(d) { return d.d1; });
    var pvals = r.digits.filter(function(d) { return d.field === 'pvalue' && d.value !== null; }).map(function(d) { return d.value; });

    // MAD for this review (may have very few digits)
    var d1Counts = new Array(9).fill(0);
    d1s.forEach(function(d) { d1Counts[d - 1]++; });
    var n = d1s.length;
    var mad = n > 0 ? computeMAD(d1Counts, bfFirst, n) : null;

    // Per-digit exact binomial p-values
    var digitPvals = [];
    for (var d = 0; d < 9; d++) {
      if (n > 0) {
        var observed = d1Counts[d];
        var expected = bfFirst[d];
        // Two-sided binomial test
        var pBinom = binomTestTwoSided(n, observed, expected);
        digitPvals.push(pBinom);
      }
    }

    // Fisher's combined p
    var fisher = fisherCombine(digitPvals);

    // Caliper contribution
    var caliperScore = 0;
    pvals.forEach(function(p) {
      if (p >= 0.04 && p < 0.05) caliperScore += 20;
      else if (p >= 0.03 && p < 0.04) caliperScore += 10;
    });

    // Suspicion score (0-100)
    var madScore = mad !== null ? Math.min(40, mad / 0.025 * 40) : 0;
    var fisherScore = fisher.p > 0 ? Math.min(40, -Math.log10(fisher.p) / 4 * 40) : 40;
    var score = Math.round(madScore + fisherScore + Math.min(20, caliperScore));

    var classification = score < 25 ? 'Conforming' :
      score < 50 ? 'Minor Deviation' :
      score < 75 ? 'Suspicious' : 'Highly Suspicious';

    return {
      review_id: r.review_id,
      analysis_name: r.analysis_name,
      k: r.k,
      n_digits: n,
      mad: mad,
      fisher_p: fisher.p,
      score: score,
      classification: classification,
      d1Counts: d1Counts,
      digits: r.digits,
    };
  });

  // BH-FDR correction
  var sorted = results.slice().sort(function(a, b) { return a.fisher_p - b.fisher_p; });
  var m = sorted.length;
  sorted.forEach(function(r, i) {
    r.fdr_p = Math.min(1, r.fisher_p * m / (i + 1));
  });
  // Enforce monotonicity
  for (var i = m - 2; i >= 0; i--) {
    sorted[i].fdr_p = Math.min(sorted[i].fdr_p, sorted[i + 1].fdr_p);
  }

  return results;
}

function binomTestTwoSided(n, k, p0) {
  // P(X >= k or X <= mirror) under Binom(n, p0)
  var pExact = 0;
  var pObs = binomPMF(n, k, p0);
  for (var x = 0; x <= n; x++) {
    var px = binomPMF(n, x, p0);
    if (px <= pObs + 1e-12) pExact += px;
  }
  return Math.min(1, pExact);
}
```

- [ ] **Step 2: Implement sortable table + accordion detail**

Table columns: Review ID, Analysis, k, Digits, MAD, Fisher p, FDR p, Score, Classification. Sortable by click. Classification badge colored. Click row to expand digit chart for that review (9 bars vs Benford expected).

- [ ] **Step 3: Verify — table shows 403 rows, sort works, accordion expands per-review chart**

- [ ] **Step 4: Commit**

```bash
cd /c/BenfordMA && git add benford-ma.html && git commit -m "feat: per-review screening table with BH-FDR + accordion detail"
```

---

### Task 7: Custom mode + export + final polish

**Files:**
- Modify: `C:\BenfordMA\benford-ma.html`

**Context:** Custom mode lets user paste CSV. Export buttons (CSV of screening table, PNG of first-digit chart). Final polish: caveat text, responsive, print styles.

- [ ] **Step 1: Implement Custom mode**

When "Custom" radio is selected:
- Show CSV paste area (textarea + Import button)
- On import: parse CSV (Study, Effect, SE, N, P-value), extract digits, run same analysis
- Hide Section 6 (per-review table) since custom is a single analysis
- Sections 1-5 render on the pasted data

```javascript
function importCustomCSV() {
  var text = document.getElementById('customCSV').value.trim();
  var lines = text.split(/\r?\n/);
  var customDigits = [];
  var customPvalues = [];

  for (var i = 1; i < lines.length; i++) { // skip header
    var cols = lines[i].split(/[,\t]/);
    if (cols.length < 5) continue;
    var fields = {
      effect: safe_float(cols[1]),
      se: safe_float(cols[2]),
      n: safe_float(cols[3]),
      pvalue: safe_float(cols[4]),
    };
    for (var fname in fields) {
      var val = fields[fname];
      if (val !== null && val !== 0) {
        var d1 = firstDigit(val), d2 = secondDigit(val), m = mantissaVal(val);
        if (d1 !== null) {
          customDigits.push({ field: fname, value: val, d1: d1, d2: d2, mantissa: m });
          if (fname === 'pvalue') customPvalues.push(val);
        }
      }
    }
  }
  // Run analysis on customDigits, render Sections 1-5
}
```

- [ ] **Step 2: Implement CSV + PNG export**

CSV: Export the per-review screening table (filtered rows).
PNG: Export the first-digit distribution chart via SVG-to-Canvas.
Same Blob URL pattern as RMST Meta / EvidenceQuality.

- [ ] **Step 3: Add specification sensitivity analysis section**

Below the main corpus analysis, add a collapsible section showing the 394K specification-level digit counts vs Benford. Uses the pre-computed `DATA.spec_sensitivity.d1_counts` and `d2_counts`. Simple bar chart + MAD.

- [ ] **Step 4: Verify full integration**

Corpus mode:
1. Summary cards correct
2. First/second digit charts with per-field breakdown
3. Mantissa arc roughly circular
4. Caliper test shows p-value clustering
5. Per-review table sortable, accordion works
6. Spec sensitivity section visible

Custom mode:
1. Paste 5 rows of CSV
2. Charts update
3. Summary cards update

Exports:
1. CSV downloads screening table
2. PNG downloads first-digit chart

Dark mode, div balance, no `</script>` in script block.

- [ ] **Step 5: Commit**

```bash
cd /c/BenfordMA && git add benford-ma.html && git commit -m "feat: custom mode, export, spec sensitivity, final polish"
```
