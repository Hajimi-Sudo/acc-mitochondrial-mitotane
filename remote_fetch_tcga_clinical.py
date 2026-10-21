import json
import pathlib
import requests

root = pathlib.Path('/root/autodl-tmp/acc-mitochondrial-mitotane')
out = root / 'data/raw/clinical'
out.mkdir(parents=True, exist_ok=True)

session = requests.Session()
session.trust_env = False
fields = ','.join([
    'case_id', 'submitter_id', 'project.project_id', 'demographic',
    'diagnoses', 'samples.submitter_id', 'samples.sample_type',
])
payload = {
    'filters': {'op': '=', 'content': {'field': 'project.project_id', 'value': 'TCGA-ACC'}},
    'fields': fields,
    'format': 'JSON',
    'size': '1000',
}
response = session.post(
    'https://api.gdc.cancer.gov/cases',
    json=payload,
    timeout=(10, 60),
)
response.raise_for_status()
data = response.json()
(out / 'TCGA_ACC_cases_clinical.json').write_text(json.dumps(data, indent=2))

hits = data.get('data', {}).get('hits', [])
flat = []
for case in hits:
    diagnosis = (case.get('diagnoses') or [{}])[0]
    demographic = case.get('demographic') or {}
    sample = (case.get('samples') or [{}])[0]
    flat.append({
        'case_id': case.get('case_id'),
        'submitter_id': case.get('submitter_id'),
        'sample_id': sample.get('submitter_id'),
        'sample_type': sample.get('sample_type'),
        'age_at_diagnosis': diagnosis.get('age_at_diagnosis'),
        'days_to_death': diagnosis.get('days_to_death'),
        'days_to_last_follow_up': diagnosis.get('days_to_last_follow_up'),
        'vital_status': diagnosis.get('vital_status'),
        'tumor_stage': diagnosis.get('tumor_stage'),
        'gender': demographic.get('gender'),
        'race': demographic.get('race'),
        'ethnicity': demographic.get('ethnicity'),
    })

with (out / 'TCGA_ACC_clinical_flat.tsv').open('w') as handle:
    handle.write('\t'.join(flat[0].keys()) + '\n')
    for row in flat:
        handle.write('\t'.join('' if v is None else str(v) for v in row.values()) + '\n')

keywords = []
def scan(value, path=''):
    if isinstance(value, dict):
        for key, item in value.items():
            scan(item, f'{path}.{key}')
    elif isinstance(value, list):
        for i, item in enumerate(value):
            scan(item, f'{path}[{i}]')
    elif any(word in path.lower() for word in ['hormone', 'functional', 'secretion', 'cortisol', 'androgen', 'mitotane', 'steroid']):
        keywords.append((path, value))

scan(data)
(out / 'TCGA_ACC_endocrine_keyword_audit.tsv').write_text(
    '\n'.join(f'{path}\t{value}' for path, value in keywords) + ('\n' if keywords else '')
)
summary = {
    'n_cases': len(hits),
    'n_flat_rows': len(flat),
    'n_endocrine_keyword_hits': len(keywords),
}
(out / 'TCGA_ACC_clinical_audit.json').write_text(json.dumps(summary, indent=2))
print(json.dumps(summary))
