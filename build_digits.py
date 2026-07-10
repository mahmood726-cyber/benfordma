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

PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_PROJECTS_ROOT = PROJECT_ROOT.parent


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
        # Prefer PredictionGap theta; fall back to FragilityAtlas median only
        # when theta is missing. Use an explicit None check, not `or`, so a
        # legitimate theta of exactly 0.0 is not silently replaced by the
        # median (a falsy 0.0 would otherwise trigger the fallback).
        theta = safe_float(pg.get('theta'))
        effect = theta if theta is not None else safe_float(fa.get('median_theta'))
        # Fields to extract from both sources
        fields = {
            'effect': effect,
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

        kval = safe_float(fa.get('k'))
        reviews.append({
            'review_id': rid,
            'analysis_name': fa.get('analysis_name', ''),
            'k': int(kval) if kval is not None else 0,
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


def resolve_paths(project_root=None, projects_root=None):
    project_root = Path(project_root).resolve() if project_root else PROJECT_ROOT
    projects_root = Path(projects_root).resolve() if projects_root else project_root.parent
    return {
        'fragility_results': projects_root / 'FragilityAtlas' / 'data' / 'output' / 'fragility_atlas_results.csv',
        'prediction_results': projects_root / 'PredictionGap' / 'data' / 'output' / 'prediction_gap_results.csv',
        'fragility_specs': projects_root / 'FragilityAtlas' / 'data' / 'output' / 'fragility_atlas_specifications.csv',
        'validation_inputs': projects_root / 'FragilityAtlas' / 'data' / 'output' / 'r_validation_inputs.json',
        'output': project_root / 'data' / 'corpus_digits.json',
    }


def main(project_root=None, projects_root=None):
    print("Building BenfordMA digit corpus...\n")
    paths = resolve_paths(project_root=project_root, projects_root=projects_root)

    # Review-level digits (primary)
    print("1. Review-level data:")
    fa = load_csv(paths['fragility_results'])
    pg = load_csv(paths['prediction_results'])
    reviews = extract_review_digits(fa, pg)

    total_digits = sum(len(r['digits']) for r in reviews)
    print(f"  -> {len(reviews)} reviews, {total_digits} digit observations\n")

    # Specification-level sensitivity (secondary)
    print("2. Specification-level sensitivity analysis:")
    spec_summary = extract_spec_digits(paths['fragility_specs'])
    print(f"  -> {spec_summary['n_values']} specification-level digit observations\n")

    # Per-study ground truth (tertiary)
    print("3. Per-study ground truth:")
    study_digits = extract_study_digits(paths['validation_inputs'])

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

    out_path = paths['output']
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
    return out_path


if __name__ == '__main__':
    main()
