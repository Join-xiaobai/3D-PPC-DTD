#!/usr/bin/env python3
"""
Step 3: Identify differentially expressed genes in cardiomyocytes 
from right ventricle of PAH vs control (GSE240921).

Input:
    data/processed/GSE240921_series_matrix.txt

Output:
    results/step3_outputs/
        ├── rv_cardiomyocyte_deg.csv
        ├── up_in_pah_rv.csv
        └── down_in_pah_rv.csv
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from io import StringIO

INPUT_FILE = Path("data/processed/GSE240921_series_matrix.txt")
OUTPUT_DIR = Path("results/step3_outputs")
OUTPUT_DIR.mkdir(exist_ok=True)


def parse_gse240921(filepath):
    """Parse GSE240921 series matrix to get cardiomyocyte expression by group."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Extract sample metadata: disease state and cell type
    samp_meta = {}
    for line in lines:
        if line.startswith('!Sample_title'):
            titles = line.strip().split('\t')[1:]
        elif line.startswith('!Sample_characteristics_ch1'):
            ch1_vals = line.strip().split('\t')[1:]
            for title, ch1 in zip(titles, ch1_vals):
                parts = ch1.split(': ')
                if len(parts) == 2:
                    key, val = parts
                    if title not in samp_meta:
                        samp_meta[title] = {}
                    samp_meta[title][key] = val
    
    # Find expression table start
    data_start = 0
    for i, line in enumerate(lines):
        if line.startswith('!series_matrix_table_begin'):
            data_start = i + 1
            break
    
    expr_lines = lines[data_start:-1]
    df = pd.read_csv(StringIO(''.join(expr_lines)), sep='\t', index_col=0)
    
    # Filter: only cardiomyocyte samples
    cardio_samples = []
    group_labels = []
    for col in df.columns:
        meta = samp_meta.get(col, {})
        cell_type = meta.get('cell type', '').lower()
        disease = meta.get('disease state', '').lower()
        
        if 'cardiomyocyte' in cell_type or 'myocyte' in cell_type:
            cardio_samples.append(col)
            if 'pulmonary arterial hypertension' in disease:
                group_labels.append('PAH')
            elif 'control' in disease:
                group_labels.append('Control')
            else:
                group_labels.append('Other')
    
    cardio_df = df[cardio_samples].copy()
    cardio_df.columns = group_labels
    
    return cardio_df


def main():
    print("Loading GSE240921 (right ventricle scRNA-seq aggregated data)...")
    expr_df = parse_gse240921(INPUT_FILE)
    
    pah_cols = [col for col in expr_df.columns if col == 'PAH']
    ctrl_cols = [col for col in expr_df.columns if col == 'Control']
    
    if len(pah_cols) < 2 or len(ctrl_cols) < 2:
        raise ValueError(f"Insufficient PAH ({len(pah_cols)}) or Control ({len(ctrl_cols)}) cardiomyocyte samples.")
    
    pah_mean = expr_df[pah_cols].mean(axis=1)
    ctrl_mean = expr_df[ctrl_cols].mean(axis=1)
    
    log2fc = np.log2((pah_mean + 1) / (ctrl_mean + 1))
    mean_expr = (pah_mean + ctrl_mean) / 2.0
    
    result = pd.DataFrame({
        'gene_symbol': log2fc.index,
        'pah_rv_cardio_mean': pah_mean.values,
        'control_rv_cardio_mean': ctrl_mean.values,
        'mean_expression': mean_expr.values,
        'log2fc_pah_vs_control_rv': log2fc.values
    })
    
    # Filter: |log2FC| > 0.58 (≈1.5x) and mean_expression > 1 (TPM-like)
    mask = (np.abs(result['log2fc_pah_vs_control_rv']) > 0.58) & (result['mean_expression'] > 1.0)
    degs = result[mask].copy()
    degs['regulation'] = degs['log2fc_pah_vs_control_rv'].apply(
        lambda x: 'up_in_pah_rv' if x > 0 else 'down_in_pah_rv'
    )
    
    # Save
    degs.to_csv(OUTPUT_DIR / "rv_cardiomyocyte_deg.csv", index=False)
    degs[degs['regulation'] == 'up_in_pah_rv'].to_csv(OUTPUT_DIR / "up_in_pah_rv.csv", index=False)
    degs[degs['regulation'] == 'down_in_pah_rv'].to_csv(OUTPUT_DIR / "down_in_pah_rv.csv", index=False)
    
    # Plot
    plt.figure(figsize=(8, 6))
    plt.scatter(result['log2fc_pah_vs_control_rv'], result['mean_expression'], s=5, alpha=0.7)
    plt.xlabel('Log2 Fold Change (PAH RV vs Control RV)')
    plt.ylabel('Mean Expression (Cardiomyocytes)')
    plt.title('Differential Expression in Right Ventricle Cardiomyocytes (GSE240921)')
    plt.axvline(-0.58, color='k', linestyle='--', linewidth=0.5)
    plt.axvline(0.58, color='k', linestyle='--', linewidth=0.5)
    plt.savefig(OUTPUT_DIR / "rv_cardiomyocyte_deg_plot.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Step 3 complete. Results saved to {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()