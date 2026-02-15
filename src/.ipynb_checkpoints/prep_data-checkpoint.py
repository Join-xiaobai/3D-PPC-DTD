#!/usr/bin/env python3
"""
prep_data.py (Enhanced Local Mode)
Preprocess locally available raw data into standardized intermediate files.

Assumes the following files exist in data/raw/:
  - GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm.gct.gz
  - GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt
  - chembl_36_activities.csv.gz
  - chembl_36_molecule_dictionary.csv.gz
  - chembl_36_target_dictionary.csv.gz
  - chembl_36_target_components.csv.gz   ← NEW!
  - uniprot_human.tsv.gz
  - GSE117261_series_matrix.txt.gz
  - GSE240921_series_matrix.txt.gz
  - GSE210248_series_matrix.txt.gz

Outputs cleaned, analysis-ready CSVs to data/processed/
"""

import gzip
import shutil
from pathlib import Path
import pandas as pd
import numpy as np

RAW_DIR = Path("data/raw")
PROC_DIR = Path("data/processed")
PROC_DIR.mkdir(parents=True, exist_ok=True)


def prepare_gtex():
    print("Processing GTEx...")
    gct_file = RAW_DIR / "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm.gct.gz"
    meta_file = RAW_DIR / "GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt"
    
    meta = pd.read_csv(meta_file, sep='\t')
    sample_to_tissue = meta.set_index('SAMPID')['SMTSD'].to_dict()
    
    with gzip.open(gct_file, 'rt') as f:
        next(f); next(f)
        expr = pd.read_csv(f, sep='\t', index_col=0)
    
    expr.index = expr.index.str.split('|').str[1]
    expr.index.name = 'gene_symbol'
    
    tissue_labels = []
    for col in expr.columns:
        samp_id = col.split('.')[0]
        tissue_labels.append(sample_to_tissue.get(samp_id, 'Other'))
    expr.columns = tissue_labels
    
    lung_tpm = expr.loc[:, expr.columns == 'Lung'].median(axis=1)
    heart_tpm = expr.loc[:, expr.columns == 'Heart - Left Ventricle'].median(axis=1)
    
    result = pd.DataFrame({
        'gene_symbol': lung_tpm.index,
        'lung_tpm': lung_tpm.values,
        'heart_tpm': heart_tpm.values
    })
    result.to_csv(PROC_DIR / "gtex_lung_heart_tpm.csv", index=False)


def prepare_uniprot():
    print("Processing UniProt...")
    uniprot_file = RAW_DIR / "uniprot_human.tsv.gz"
    df = pd.read_csv(uniprot_file, sep='\t')
    df = df.rename(columns={'Entry': 'uniprot_id', 'Gene Names': 'gene_symbol'})
    df['gene_symbol'] = df['gene_symbol'].astype(str).str.split().str[0]
    mapping = df[['uniprot_id', 'gene_symbol']].dropna().drop_duplicates()
    mapping.to_csv(PROC_DIR / "uniprot_human_mapping.csv", index=False)


def prepare_chembl():
    print("Processing ChEMBL with gene symbol mapping...")
    # Load tables
    activities = pd.read_csv(RAW_DIR / "chembl_36_activities.csv.gz")
    molecules = pd.read_csv(RAW_DIR / "chembl_36_molecule_dictionary.csv.gz")
    targets = pd.read_csv(RAW_DIR / "chembl_36_target_dictionary.csv.gz")
    target_components = pd.read_csv(RAW_DIR / "chembl_36_target_components.csv.gz")  # NEW!
    uniprot_map = pd.read_csv(RAW_DIR / "uniprot_human.tsv.gz", sep='\t')
    uniprot_map = uniprot_map.rename(columns={'Entry': 'uniprot_id', 'Gene Names': 'gene_symbol'})
    uniprot_map['gene_symbol'] = uniprot_map['gene_symbol'].astype(str).str.split().str[0]
    uniprot_to_gene = uniprot_map.set_index('uniprot_id')['gene_symbol'].to_dict()

    # Human targets
    human_targets = targets[targets['organism'] == 'Homo sapiens'][['tid', 'chembl_id']].rename(
        columns={'chembl_id': 'target_chembl_id'}
    )

    # Map tid → UniProt → gene_symbol
    tid_to_uniprot = target_components[target_components['component_type'] == 'PROTEIN'][
        ['tid', 'component_id']
    ].rename(columns={'component_id': 'uniprot_id'})
    
    tid_to_gene = tid_to_uniprot.copy()
    tid_to_gene['target_gene_symbol'] = tid_to_gene['uniprot_id'].map(uniprot_to_gene)
    tid_to_gene = tid_to_gene.dropna(subset=['target_gene_symbol'])[['tid', 'target_gene_symbol']]

    # Approved drugs
    approved_mols = molecules[molecules['max_phase'] >= 4][['chembl_id', 'pref_name']]

    # High-confidence activities
    high_conf_act = activities[
        (activities['pchembl_value'] >= 6) &
        (activities['standard_type'].isin(['IC50', 'EC50', 'Ki', 'Kd'])) &
        (activities['relation'] == '=')
    ][['molecule_chembl_id', 'tid', 'pchembl_value']]

    # Merge all
    df = (
        high_conf_act
        .merge(human_targets, on='tid', how='inner')
        .merge(tid_to_gene, on='tid', how='inner')
        .merge(approved_mols, left_on='molecule_chembl_id', right_on='chembl_id', how='inner')
    )

    final = df[[
        'molecule_chembl_id', 'pref_name',
        'target_chembl_id', 'target_gene_symbol',
        'pchembl_value'
    ]].drop_duplicates().rename(columns={'pref_name': 'drug_name'})

    final.to_csv(PROC_DIR / "chembl_drug_target_pairs.csv", index=False)


def extract_geo_matrices():
    print("Extracting GEO matrices...")
    for gse in ["GSE117261", "GSE240921", "GSE210248"]:
        gz_file = RAW_DIR / f"{gse}_series_matrix.txt.gz"
        txt_file = PROC_DIR / f"{gse}_series_matrix.txt"
        if not txt_file.exists():
            with gzip.open(gz_file, 'rt') as f_in, open(txt_file, 'w') as f_out:
                shutil.copyfileobj(f_in, f_out)


def main():
    print("Starting enhanced local data preprocessing...\n")
    
    prepare_gtex()
    prepare_uniprot()
    prepare_chembl()          # Now includes target_gene_symbol
    extract_geo_matrices()
    
    print("\n✅ All preprocessing complete.")
    print(f"Intermediate files saved to: {PROC_DIR.resolve()}")


if __name__ == "__main__":
    main()