Mahmood Ahmad
Tahir Heart Institute
author@example.com

Benford Screening of 1.2 Million Meta-Analytic Values Finds No Corpus-Level Digit Anomaly

Do the numerical outputs of large-scale meta-analytic workflows follow the digit patterns expected under Benford's law, or do they show overall corpus-level anomalies that might suggest data integrity concerns? We extracted first and second significant digits from 1,175,056 values across six numeric fields in 403 Cochrane review specifications from the Pairwise70 corpus. Digit frequencies were compared with Benford expectations using mean absolute deviation, chi-squared goodness-of-fit testing, and mantissa arc uniformity analysis applied to each field separately. The corpus-level first-digit mean absolute deviation was 0.013 (95% CI 0.011 to 0.015), remaining within accepted conformity thresholds despite the very large sample size. Sensitivity analyses across numeric fields and review subgroups showed consistent conformity, and no single review exceeded the critical non-conformance threshold after Bonferroni correction. At corpus level, these Cochrane outputs do not show digit anomalies suggesting systematic fabrication or major computational distortion. Benford screening remains indirect and cannot detect fabrication strategies that deliberately preserve expected digit frequencies.

Outside Notes

Type: methods
Primary estimand: Mean absolute deviation (MAD)
App: BenfordMA v1.0
Data: Pairwise70 specification-level data (403 reviews, 1.18M values)
Code: https://github.com/mahmood726-cyber/benfordma
Version: 1.0
Validation: Author reviewed draft

References

1. Carlisle JB. Data fabrication and other reasons for non-random sampling in 5087 randomised, controlled trials in anaesthetic and general medical journals. Anaesthesia. 2017;72(8):944-952.
2. Brown NJL, Heathers JAJ. The GRIM test: a simple technique detects numerous anomalies in the reporting of results in psychology. Soc Psychol Personal Sci. 2017;8(4):363-369.
3. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.

AI Disclosure

This work represents a compiler-generated evidence micro-publication (i.e., a structured, pipeline-based synthesis output). AI (Claude, Anthropic) was used as a constrained synthesis engine operating on structured inputs and predefined rules for infrastructure generation, not as an autonomous author. The 156-word body was written and verified by the author, who takes full responsibility for the content. This disclosure follows ICMJE recommendations (2023) that AI tools do not meet authorship criteria, COPE guidance on transparency in AI-assisted research, and WAME recommendations requiring disclosure of AI use. All analysis code, data, and versioned evidence capsules (TruthCert) are archived for independent verification.
