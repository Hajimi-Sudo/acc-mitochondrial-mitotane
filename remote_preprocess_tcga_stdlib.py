import csv
import gzip
import json
import pathlib

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

with (OUT / 'TCGA_ACC_file_sample_manifest.tsv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter='\t')
    writer.writeheader()
    writer.writerows(rows)

gene_order = []
gene_names = {}
sample_values = {}
for index, row in enumerate(rows, 1):
    values = {}
    with open(row['file'], encoding='utf-8') as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[:4] == ['gene_id', 'gene_name', 'gene_type', 'unstranded']:
                break
        for fields in csv.reader(handle, delimiter='\t'):
            if len(fields) < 4 or fields[2] != 'protein_coding':
                continue
            gene_id = fields[0].split('.', 1)[0]
            if gene_id not in gene_names:
                gene_order.append(gene_id)
                gene_names[gene_id] = fields[1]
            values[gene_id] = fields[3]
    sample_values[row['sample_id']] = values
    print(f'processed={index}/{len(rows)} sample={row["sample_id"]}', flush=True)

matrix_path = OUT / 'TCGA_ACC_protein_coding_unstranded_counts.tsv.gz'
with gzip.open(matrix_path, 'wt', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['gene_id', 'gene_name'] + [row['sample_id'] for row in rows])
    for gene_id in gene_order:
        writer.writerow([gene_id, gene_names[gene_id]] + [sample_values[row['sample_id']].get(gene_id, '0') for row in rows])

qc = {
    'n_files': len(rows),
    'n_samples': len(rows),
    'n_protein_coding_genes': len(gene_order),
    'sample_manifest': str(OUT / 'TCGA_ACC_file_sample_manifest.tsv'),
    'counts_matrix': str(matrix_path),
}
(OUT / 'TCGA_ACC_preprocess_qc.json').write_text(json.dumps(qc, indent=2))
print(json.dumps(qc))
