import json
import pathlib
import pandas as pd

ROOT = pathlib.Path('/root/autodl-tmp/acc-mitochondrial-mitotane')
RAW = ROOT / 'data/raw/tcga_star_counts'
OUT = ROOT / 'data/processed'
OUT.mkdir(parents=True, exist_ok=True)

inventory = json.loads(
    (ROOT / 'manifests/tcga_acc_star_counts_inventory.json').read_text()
)['data']['hits']

rows = []
for record in inventory:
    file_id = record['file_id']
    filename = pathlib.Path(record.get('file_name') or f'{file_id}.dat').name
    path = RAW / filename
    if not path.exists():
        path = RAW / f'{file_id}__{filename}'
    case = record.get('cases', [{}])[0]
    sample = case.get('samples', [{}])[0]
    rows.append({
        'file_id': file_id,
        'file': str(path),
        'case_id': case.get('submitter_id'),
        'sample_id': sample.get('submitter_id'),
        'sample_type': sample.get('sample_type'),
    })

manifest = pd.DataFrame(rows)
manifest.to_csv(
    OUT / 'TCGA_ACC_file_sample_manifest.tsv', sep='\t', index=False
)

frames = []
for row in rows:
    raw = pd.read_csv(
        row['file'],
        sep='\t',
        comment='#',
        usecols=['gene_id', 'gene_name', 'gene_type', 'unstranded'],
    )
    raw = raw[raw['gene_type'].eq('protein_coding')].copy()
    raw['gene_id'] = raw['gene_id'].str.replace(r'\..*$', '', regex=True)
    raw = raw.drop_duplicates('gene_id').set_index('gene_id')
    frame = raw[['gene_name', 'unstranded']].rename(
        columns={'unstranded': row['sample_id']}
    )
    frames.append(frame)

matrix = frames[0][['gene_name']]
for frame in frames:
    matrix = matrix.join(frame.drop(columns=['gene_name']), how='outer')
matrix = matrix.reset_index().rename(columns={'index': 'gene_id'})
matrix.to_csv(
    OUT / 'TCGA_ACC_protein_coding_unstranded_counts.tsv.gz',
    sep='\t', index=False, compression='gzip'
)

qc = {
    'n_files': len(rows),
    'n_samples': int(matrix.shape[1] - 2),
    'n_protein_coding_genes': int(matrix.shape[0]),
    'sample_manifest': str(OUT / 'TCGA_ACC_file_sample_manifest.tsv'),
    'counts_matrix': str(OUT / 'TCGA_ACC_protein_coding_unstranded_counts.tsv.gz'),
}
(OUT / 'TCGA_ACC_preprocess_qc.json').write_text(json.dumps(qc, indent=2))
print(json.dumps(qc))
