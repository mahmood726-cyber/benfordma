import csv
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from build_digits import (
    extract_review_digits,
    first_digit,
    main,
    mantissa,
    resolve_paths,
    second_digit,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_main_uses_repo_relative_sibling_projects(tmp_path):
    projects_root = tmp_path / 'projects'
    project_root = projects_root / 'BenfordMA'
    project_root.mkdir(parents=True)

    paths = resolve_paths(project_root=project_root, projects_root=projects_root)

    write_csv(
        paths['fragility_results'],
        [{
            'review_id': 'CD000001',
            'analysis_name': 'Synthetic analysis',
            'k': '4',
            'median_theta': '0.45',
            'iqr_theta': '0.12',
        }],
        ['review_id', 'analysis_name', 'k', 'median_theta', 'iqr_theta'],
    )
    write_csv(
        paths['prediction_results'],
        [{
            'review_id': 'CD000001',
            'theta': '0.55',
            'p_value': '0.03',
            'tau2': '0.01',
            'I2': '42.0',
        }],
        ['review_id', 'theta', 'p_value', 'tau2', 'I2'],
    )
    write_csv(
        paths['fragility_specs'],
        [{'theta': '0.55', 'se_theta': '0.12', 'p_value': '0.03'}],
        ['theta', 'se_theta', 'p_value'],
    )
    paths['validation_inputs'].parent.mkdir(parents=True, exist_ok=True)
    paths['validation_inputs'].write_text(
        json.dumps([{
            'review_id': 'CD000001',
            'yi': [0.1],
            'sei': [0.2],
        }]),
        encoding='utf-8',
    )

    out_path = main(project_root=project_root, projects_root=projects_root)

    payload = json.loads(out_path.read_text(encoding='utf-8'))
    assert out_path == project_root / 'data' / 'corpus_digits.json'
    assert payload['corpus_summary'] == {
        'n_reviews': 1,
        'n_review_digits': 6,
        'n_spec_digits': 3,
        'n_study_digits': 2,
        'fields': ['effect', 'se', 'k', 'pvalue', 'tau2', 'I2'],
    }
    assert payload['reviews'][0]['review_id'] == 'CD000001'
    assert payload['reviews'][0]['digits'][0]['field'] == 'effect'


# ---- build_digits digit-extraction edge cases (F3 coverage) ----

def test_first_second_digit_and_mantissa_edge_cases():
    # None / 0 / non-finite -> None
    for bad in (None, 0, float('inf'), float('nan')):
        assert first_digit(bad) is None
        assert second_digit(bad) is None
        assert mantissa(bad) is None

    # Values below 1 are still handled and the mantissa stays in [0, 1).
    for val in (0.45, 0.12, 4.5, 45.0, -0.045, 1.0):
        m = mantissa(val)
        assert m is not None
        assert 0.0 <= m < 1.0, (val, m)

    # First/second digit are sign- and scale-invariant.
    assert first_digit(0.45) == 4
    assert second_digit(0.45) == 5
    assert first_digit(-45.0) == 4
    assert second_digit(-45.0) == 5

    # Mantissa matches the closed form for a sub-1 value (the F1 regression:
    # log10(0.45) % 1 must be the NON-negative fractional part ~0.6532).
    assert math.isclose(mantissa(0.45), 0.6532125137753437, rel_tol=1e-9)


def test_extract_review_digits_tolerates_blank_k(tmp_path):
    # F2 regression: a present-but-blank k cell must not crash the build.
    fa_rows = [{
        'review_id': 'CD000002',
        'analysis_name': 'Blank k row',
        'k': '',            # blank cell -> float('') would raise ValueError
        'median_theta': '0.45',
        'iqr_theta': '0.12',
    }]
    pg_rows = [{
        'review_id': 'CD000002',
        'theta': '0.55',
        'p_value': '0.03',
        'tau2': '0.01',
        'I2': '42.0',
    }]

    reviews = extract_review_digits(fa_rows, pg_rows)

    assert len(reviews) == 1
    assert reviews[0]['review_id'] == 'CD000002'
    assert reviews[0]['k'] == 0  # blank -> safe default, no crash


# ---- JS statistical core: mantissaVal sign normalization (F1 regression) ----

def test_js_mantissa_val_normalized_to_unit_interval():
    node = shutil.which('node')
    html = (REPO_ROOT / 'benford-ma.html').read_text(encoding='utf-8')

    # Source-level guard so the fix is locked even when Node is unavailable.
    assert 'm < 0 ? m + 1 : m' in html, 'mantissaVal must normalize negatives into [0,1)'

    if node is None:
        return  # Node not installed in this environment; source check above suffices.

    # Extract the mantissaVal function body and exercise it under Node.
    start = html.index('function mantissaVal(')
    end = html.index('\n}', start) + 2
    fn_src = html[start:end]
    harness = (
        fn_src
        + "\nvar cases=[[0.45,0.6532125137753437],[0.12,0.07918124604762478],"
          "[4.5,0.6532125137753437]];\n"
        + "for (var i=0;i<cases.length;i++){var got=mantissaVal(cases[i][0]);"
          "if(got<0||got>=1){console.error('OUT_OF_RANGE',cases[i][0],got);process.exit(2);}"
          "if(Math.abs(got-cases[i][1])>1e-9){console.error('MISMATCH',cases[i][0],got);process.exit(3);}}\n"
        + "console.log('OK');\n"
    )
    result = subprocess.run(
        [node, '-e', harness], capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert 'OK' in result.stdout
