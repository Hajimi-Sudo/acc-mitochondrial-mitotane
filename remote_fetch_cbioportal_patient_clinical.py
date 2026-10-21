import concurrent.futures
import json
import pathlib
import requests

root = pathlib.Path('/root/autodl-tmp/acc-mitochondrial-mitotane')
out = root / 'data/raw/clinical'
out.mkdir(parents=True, exist_ok=True)
patients = requests.get(
    'https://www.cbioportal.org/api/studies/acc_tcga/patients',
    params={'projection': 'DETAILED', 'pageSize': 1000},
    timeout=(8, 60),
).json()

def fetch(patient):
    pid = patient['patientId']
    url = f'https://www.cbioportal.org/api/studies/acc_tcga/patients/{pid}/clinical-data'
    for attempt in range(4):
        try:
            session = requests.Session()
            session.trust_env = False
            response = session.get(url, timeout=(8, 30))
            response.raise_for_status()
            return pid, response.json(), ''
        except Exception as exc:
            if attempt == 3:
                return pid, [], repr(exc)
    return pid, [], 'unreachable'

records = []
errors = []
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
    for pid, data, error in pool.map(fetch, patients):
        records.extend(data)
        if error:
            errors.append({'patientId': pid, 'error': error})
        print(f'patient={pid} records={len(data)} error={error}', flush=True)

(out / 'cbioportal_acc_tcga_patient_clinical_data.json').write_text(
    json.dumps(records, indent=2)
)
(out / 'cbioportal_acc_tcga_patient_clinical_errors.json').write_text(
    json.dumps(errors, indent=2)
)
print(json.dumps({'n_patients': len(patients), 'n_records': len(records), 'n_errors': len(errors)}))
