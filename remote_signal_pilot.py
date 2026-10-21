import json
import pathlib
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

ROOT = pathlib.Path('/root/autodl-tmp/acc-mitochondrial-mitotane')
OUT = ROOT / 'results/pilot'
OUT.mkdir(parents=True, exist_ok=True)
counts = pd.read_csv(ROOT / 'data/processed/TCGA_ACC_protein_coding_unstranded_counts.tsv.gz', sep='\t').set_index('gene_name')
counts = counts.drop(columns=['gene_id'])
counts = counts.loc[~counts.index.duplicated()]
expr = np.log2(counts.astype(float) + 1.0)
sets = {
    'OXPHOS': pathlib.Path(ROOT / 'data/raw/gene_sets/hallmark_oxphos.dat').read_text().splitlines()[2:],
    'ROS': pathlib.Path(ROOT / 'data/raw/gene_sets/hallmark_ros.dat').read_text().splitlines()[2:],
    'G2M': pathlib.Path(ROOT / 'data/raw/gene_sets/hallmark_g2m.dat').read_text().splitlines()[2:],
}
scores = pd.DataFrame(index=expr.columns)
audit = {}
for name, genes in sets.items():
    genes = sorted(set(x.strip().upper() for x in genes if x.strip()) & set(expr.index))
    audit[name] = len(genes)
    z = expr.loc[genes].sub(expr.loc[genes].mean(axis=1), axis=0).div(expr.loc[genes].std(axis=1).replace(0, np.nan), axis=0)
    scores[name] = z.mean(axis=0)
scores['MITO_COMPOSITE'] = scores['OXPHOS']
scores['MITO_ROS'] = scores[['OXPHOS', 'ROS']].mean(axis=1)
scores['G2M_CONTROL'] = scores['G2M']
cluster_input = scores[['MITO_COMPOSITE', 'ROS']]
labels = KMeans(n_clusters=2, n_init=100, random_state=20261020).fit_predict(StandardScaler().fit_transform(cluster_input))
means = scores.assign(label=labels).groupby('label')['MITO_COMPOSITE'].mean()
high = int(means.idxmax())
scores['MITO_CLUSTER'] = np.where(labels == high, 'MITO_HIGH', 'MITO_LOW')
clinical = json.loads((ROOT / 'data/raw/clinical/cbioportal_acc_tcga_patient_clinical_data.json').read_text())
attr = {}
for record in clinical:
    attr.setdefault(record['patientId'], {})[record['clinicalAttributeId']] = record.get('value')
manifest = pd.read_csv(ROOT / 'data/processed/TCGA_ACC_file_sample_manifest.tsv', sep='\t')
rows = []
for _, row in manifest.iterrows():
    item = {'sample_id': row['sample_id'], 'patient_id': row['case_id']}
    item.update(attr.get(row['case_id'], {}))
    rows.append(item)
clinical_df = pd.DataFrame(rows).set_index('sample_id')
data = scores.join(clinical_df, how='left')
tests = {}
for col in ['HISTORY_ADRENAL_HORMONE_EXCESS', 'PHARM_TX_MITOTANE_INDICATOR', 'PHARM_TX_MITOTANE_THERAPUTIC_LEVELS']:
    x = data.dropna(subset=['MITO_COMPOSITE', col]).copy()
    high_group = x['MITO_COMPOSITE'] >= x['MITO_COMPOSITE'].median()
    table = pd.crosstab(high_group, x[col].astype(str))
    if table.shape == (2, 2):
        odds, p = fisher_exact(table.to_numpy())
        tests[col] = {'n': len(x), 'table': table.to_dict(), 'odds_ratio': float(odds), 'p': float(p)}
    else:
        tests[col] = {'n': len(x), 'table': table.to_dict()}
for col in ['WEISS_SCORE_OVERALL', 'OS_MONTHS', 'DFS_MONTHS']:
    values = pd.to_numeric(data.get(col, pd.Series(index=data.index)), errors='coerce')
    keep = values.notna()
    if keep.sum() >= 3:
        rho, p = spearmanr(data.loc[keep, 'MITO_COMPOSITE'], values[keep])
        tests[col] = {'n': int(keep.sum()), 'rho': float(rho), 'p': float(p)}
results = {'stage': 'remote_signal_pilot', 'n_samples': len(data), 'cluster_counts': scores['MITO_CLUSTER'].value_counts().to_dict(), 'gene_set_overlap': audit, 'mito_ros_rho': float(spearmanr(scores['MITO_COMPOSITE'], scores['ROS']).statistic), 'mito_g2m_rho': float(spearmanr(scores['MITO_COMPOSITE'], scores['G2M_CONTROL']).statistic), 'tests': tests, 'claim_boundary': 'association only; mitotane fields are exposure/therapeutic-level proxies, not direct response labels'}
data.to_csv(OUT / 'tcga_acc_mitochondrial_signal_scores.tsv.gz', sep='\t', compression='gzip')
(OUT / 'tcga_acc_mitochondrial_signal_pilot.json').write_text(json.dumps(results, indent=2, default=str))
print(json.dumps(results, indent=2, default=str))
