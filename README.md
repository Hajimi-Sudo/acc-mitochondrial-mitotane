# Mitochondrial oxidative-metabolism analysis code

This directory contains scripts for public-data acquisition, preprocessing,
statistical analysis, robustness assessment, and figure-data generation. They
are designed for execution on the isolated remote project root
`/root/autodl-tmp/acc-mitochondrial-mitotane`.

Raw expression matrices, patient-level clinical records, controlled-access
files, credentials, and remote result directories are not part of this code
release.

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
- `remote_audit_gsa_access.py`
- `remote_search_public_drug_sensitivity.py`
- `remote_audit_new_mitotane_paper.py`

Clinical and expression preprocessing:

- `remote_fetch_tcga_clinical.py`
- `remote_fetch_cbioportal_patient_clinical.py`
- `remote_preprocess_tcga.py`
- `remote_preprocess_tcga_stdlib.py`
- `remote_formal_mito_analysis.py`
- `remote_discovery_module.py`
- `remote_signal_pilot.py`

Public model and external-context analyses:

- `remote_gse140818_mito_analysis.py`
- `remote_gse140818_effect_stats.py`
- `remote_gse140818_adjusted_stats.py`
- `remote_gse140818_factorial_analysis.py`
- `remote_gse73417_xenograft_mito_analysis.py`
- `remote_external_gse19750.py`
- `remote_add_gse19750_package.py`

Robustness, multiplicity, and evidence packaging:

- `remote_robustness_analysis.py`
- `remote_multiple_testing_audit.py`
- `remote_sensitivity_gate.py`
- `remote_build_evidence_manifest.py`
- `remote_make_robustness_package.py`
- `remote_make_publication_package.py`

Figure rendering:

- `remote_render_jpg_figures.py`
- `render_aggregate_figures.py`

## Execution order

1. Run remote identity, storage, access, and public-data audits.
2. Acquire and preprocess data under the remote project root.
3. Run mitochondrial scoring and state-construction analysis.
4. Run external expression, endocrine, HAC15, and xenograft context analyses.
5. Run robustness and multiple-testing audits.
6. Build aggregate publication outputs and the evidence manifest.

The exact paths and cohort roles are recorded in the project contracts under
`00_contract/`. The evidence ceiling remains hypothesis generation. These
scripts do not create a patient-level mitotane-response classifier without a
sample-matched response endpoint.

## Runtime requirements

The analysis uses Python 3.10 or later with the packages listed in
`requirements.txt`. Network access is required only for remote acquisition and
metadata-audit scripts. Analysis scripts read and write only under the
configured remote project root.
