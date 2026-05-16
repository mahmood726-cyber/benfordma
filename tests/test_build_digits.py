import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from build_digits import main, resolve_paths


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
