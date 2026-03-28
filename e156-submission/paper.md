Mahmood Ahmad
Tahir Heart Institute
mahmood.ahmad2@nhs.net

Benford Screening of 1.2 Million Meta-Analytic Values Finds No Corpus-Level Digit Anomaly

Do the numerical outputs of large-scale meta-analytic workflows follow the digit patterns expected under Benford's law, or do they show corpus-level anomalies? We extracted first and second significant digits from 1,175,056 values across six numeric fields in 403 Cochrane review specifications from the Pairwise70 corpus. Digit frequencies were compared with Benford expectations using mean absolute deviation, chi-squared goodness-of-fit testing, and mantissa arc uniformity analysis. The corpus-level first-digit mean absolute deviation was 0.013 (95% CI 0.011 to 0.015), which remained within accepted conformity thresholds despite the very large sample size. Sensitivity analyses across numeric fields and review subgroups showed the same pattern, and no single review exceeded the critical non-conformance threshold after Bonferroni correction. At corpus level, these Cochrane outputs do not show the digit anomalies that would suggest systematic fabrication or major computational distortion. Benford screening remains indirect and would miss fabrication strategies that preserve expected digit frequencies.

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

This work represents a compiler-generated evidence micro-publication (i.e., a structured, pipeline-based synthesis output). AI is used as a constrained synthesis engine operating on structured inputs and predefined rules, rather than as an autonomous author. Deterministic components of the pipeline, together with versioned, reproducible evidence capsules (TruthCert), are designed to support transparent and auditable outputs. All results and text were reviewed and verified by the author, who takes full responsibility for the content. The workflow operationalises key transparency and reporting principles consistent with CONSORT-AI/SPIRIT-AI, including explicit input specification, predefined schemas, logged human-AI interaction, and reproducible outputs.
