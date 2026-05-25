# Benford Screening of 1.2 Million Meta-Analytic Values Finds No Corpus-Level Digit Anomaly

[![ci](https://github.com/mahmood726-cyber/benfordma/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/mahmood726-cyber/benfordma/actions/workflows/ci.yml) [![codeql](https://github.com/mahmood726-cyber/benfordma/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/mahmood726-cyber/benfordma/actions/workflows/codeql.yml) [![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE) [![python: 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)

Do the numerical outputs of large-scale meta-analytic workflows follow the digit patterns expected under Benford's law, or do they show overall corpus-level anomalies that might suggest data integrity concerns? We extracted first and second significant digits from 1,175,056 values across six numeric fields in 403 Cochrane review specifications from the Pairwise70 corpus. Digit frequencies were compared with Benford expectations using mean absolute deviation, chi-squared goodness-of-fit testing, and mantissa arc uniformity analysis applied to each field separately. The corpus-level first-digit mean absolute deviation was 0.013 (95% CI 0.011 to 0.015), remaining within accepted conformity thresholds despite the very large sample size. Sensitivity analyses across numeric fields and review subgroups showed consistent conformity, and no single review exceeded the critical non-conformance threshold after Bonferroni correction. At corpus level, these Cochrane outputs do not show digit anomalies suggesting systematic fabrication or major computational distortion. Benford screening remains indirect and cannot detect fabrication strategies that deliberately preserve expected digit frequencies.

**Live dashboard:** <https://mahmood726-cyber.github.io/benfordma/>

## Run

Open `benford-ma.html` (or `index.html`) in any modern browser. No build step.

For local development:

```bash
python -m http.server 8000
# then open http://localhost:8000/
```

## Test

```bash
python -m pytest -q
```

The suite under `tests/` includes 2 test file(s).

## Repo layout

| Path | Purpose |
|---|---|
| `benford-ma.html` | the dashboard (main artifact) |
| `index.html` | landing page |
| `tests/` | pytest tests |
| `e156-submission/` | E156 micro-paper bundle |
| `E156-PROTOCOL.md` | project metadata (E156 entry #11) |

## License

See `LICENSE` (MIT).
