# ACC mitochondrial oxidative-metabolism analysis

This repository contains the remote-only scripts for scoring mitochondrial
oxidative metabolism and analyzing public adrenocortical carcinoma expression
data in treatment-related model contexts.

The repository contains code only. It does not contain raw expression matrices,
patient-level clinical records, controlled-access files, credentials, private
execution outputs, or manuscript files.

## Execution boundary

All data acquisition, preprocessing, statistical analysis, and figure
generation must run under an isolated remote project root. The scripts use the
default path `/root/autodl-tmp/acc-mitochondrial-mitotane`, which can be
overridden only where a script explicitly exposes a project-root option or
environment variable.

Do not run the data acquisition or analysis scripts on a local machine that is
not approved for the data. Do not place restricted or patient-level data in a
Git repository.

## Data sources

The scripts use public resources identified by accession: TCGA-ACC, GSE76019,
GSE140818, GSE186870, GSE73417, and GSE19750. HRA006596 is controlled access
and is not included or downloaded by this repository.

## Script groups

Data access and audit:

- `remote_pilot_audit.sh`
- `remote_data_availability_audit.sh`
- `remote_public_data_availability_audit.py`
- `remote_public_geo_qc.py`
- `remote_download_inventory.py`
- `remote_download_aria2.sh`
- `make_aria2_list.py`
- `remote_validate_tcga_download.py`
- `remote_search_public_drug_sensitivity.py`

Preprocessing and discovery analysis:

- `remote_fetch_tcga_clinical.py`
- `remote_fetch_cbioportal_patient_clinical.py`
- `remote_preprocess_tcga.py`
- `remote_preprocess_tcga_stdlib.py`
- `remote_formal_mito_analysis.py`
- `remote_discovery_module.py`
- `remote_signal_pilot.py`

External and model-system analyses:

- `remote_external_gse19750.py`
- `remote_gse140818_mito_analysis.py`
- `remote_gse140818_effect_stats.py`
- `remote_gse140818_adjusted_stats.py`
- `remote_gse140818_factorial_analysis.py`
- `remote_gse73417_xenograft_mito_analysis.py`

Robustness and figure generation:

- `remote_robustness_analysis.py`
- `remote_multiple_testing_audit.py`
- `remote_sensitivity_gate.py`
- `remote_render_jpg_figures.py`

## Suggested execution order

1. Run the remote identity, storage, access, and public-data audits.
2. Acquire and preprocess public data under the isolated remote root.
3. Run mitochondrial scoring and state construction.
4. Run external expression, endocrine, HAC15, and xenograft context analyses.
5. Run robustness and multiple-testing audits.
6. Render aggregate figures from recorded analysis outputs.

The scripts support hypothesis generation and model-system context. They do not
establish a patient-level mitotane-response classifier without sample-matched
pretreatment molecular data and a direct response endpoint.

## Requirements

Python 3.10 or later is recommended. Install the packages in
`requirements.txt` in the remote environment. Network access is required only
for remote metadata and public-data acquisition scripts.

## License

The code is released under the MIT License. See `LICENSE`.
