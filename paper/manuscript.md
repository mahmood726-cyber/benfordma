# Benford Screening of 1.2 Million Meta-Analytic Values Finds No Corpus-Level Digit Anomaly

## Overview

First-digit distributions of 1.2 million meta-analytic values conform to Benford expectations with MAD = 0.013. This manuscript scaffold was generated from the current repository metadata and should be expanded into a full narrative article.

## Study Profile

Type: methods
Primary estimand: Mean absolute deviation (MAD)
App: BenfordMA v1.0
Data: Pairwise70 specification-level data (403 reviews, 1.18M values)
Code: https://github.com/mahmood726-cyber/benfordma

## E156 Capsule

Do the numerical outputs of large-scale meta-analytic workflows follow the digit patterns expected under Benford's law, or do they show overall corpus-level anomalies that might suggest data integrity concerns? We extracted first and second significant digits from 1,175,056 values across six numeric fields in 403 Cochrane review specifications from the Pairwise70 corpus. Digit frequencies were compared with Benford expectations using mean absolute deviation, chi-squared goodness-of-fit testing, and mantissa arc uniformity analysis applied to each field separately. The corpus-level first-digit mean absolute deviation was 0.013 (95% CI 0.011 to 0.015), remaining within accepted conformity thresholds despite the very large sample size. Sensitivity analyses across numeric fields and review subgroups showed consistent conformity, and no single review exceeded the critical non-conformance threshold after Bonferroni correction. At corpus level, these Cochrane outputs do not show digit anomalies suggesting systematic fabrication or major computational distortion. Benford screening remains indirect and cannot detect fabrication strategies that deliberately preserve expected digit frequencies.

## Expansion Targets

1. Expand the background and rationale into a full introduction.
2. Translate the E156 capsule into detailed methods, results, and discussion sections.
3. Add figures, tables, and a submission-ready reference narrative around the existing evidence object.
