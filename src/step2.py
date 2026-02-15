#!/usr/bin/env python3
"""
Step 2: Identify differentially expressed genes (DEGs) in PAH lung tissue 
using GEO dataset GSE117261.

Input:
    data/processed/GSE117261_series_matrix.txt  (from prep_data.py)

Output:
    results/step2_outputs/
        ├── pah_lung_deg.csv
        ├── upregulated_in_pah.csv
        └── downregulated_in_pah.csv
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Paths
INPUT_FILE = Path("data/processed/GSE117261_series_matrix.txt")
OUTPUT_DIR = Path("results/step2_outputs")
OUTPUT_DIR.mkdir(exist_ok=True)


def parse_geo_series_matrix(filepath):
    """Parse GEO series matrix to expression DataFrame and sample metadata."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find data start line
    data_start = 0
    for i, line in enumerate(lines):
        if line.startswith('!series_matrix_table_begin'):
            data_start = i + 1
            break
    
    # Read expression data
    expr_lines = lines[data_start:-1]  # exclude last end line
    from io import StringIO
    df = pd.read_csv(StringIO(''.join(expr_lines)), sep='\t', index_col=0)
    
    # Extract sample characteristics (disease state)
    sample_groups = {}
    for line in lines:
        if line.startswith('!Sample_characteristics_ch1'):
            parts = line.strip().split('\t')
            key_val = parts[1].split(': ')
            if len(key_val) == 2:
                key, val = key_val
                if key == 'disease state':
                    samp_id = parts[0].split('_')[-1]
                    sample_groups[samp_id] = val
    
    # Map columns to disease state
    col_labels = []
    for col in df.columns:
        state = sample_groups.get(col, 'Unknown')
        col_labels.append(state)
    
    df.columns = col_labels
    return df


def main():
    print("Loading GSE117261 data...")
    expr_df = parse_geo_series_matrix(INPUT_FILE)
    
    # Separate PAH and Control
    pah_cols = [col for col in expr_df.columns if 'pulmonary arterial hypertension' in str(col).lower()]
    ctrl_cols = [col for col in expr_df.columns if 'control' in str(col).lower()]
    
    if not pah_cols or not ctrl_cols:
        raise ValueError("Could not identify PAH or control samples in GSE117261.")
    
    pah_expr = expr_df[pah_cols].mean(axis=1)
    ctrl_expr = expr_df[ctrl_cols].mean(axis=1)
    
    # Compute log2 fold-change (PAH / Control)
    log2fc = np.log2((pah_expr + 1) / (ctrl_expr + 1))
    mean_expr = (pah_expr + ctrl_expr) / 2.0
    
    result = pd.DataFrame({
        'gene_symbol': log2fc.index,
        'pah_mean': pah_expr.values,
        'control_mean': ctrl_expr.values,
        'mean_expression': mean_expr.values,
        'log2fc_pah_vs_control': log2fc.values
    })
    
    # Filter significant DEGs: |log2FC| > 1 and mean_expression > 5
    mask = (np.abs(result['log2fc_pah_vs_control']) > 1.0) & (result['mean_expression'] > 5.0)
    degs = result[mask].copy()
    degs['regulation'] = degs['log2fc_pah_vs_control'].apply(
        lambda x: 'up_in_pah' if x > 0 else 'down_in_pah'
    )
    
    # Save results
    degs.to_csv(OUTPUT_DIR / "pah_lung_deg.csv", index=False)
    
    up_genes = degs[degs['regulation'] == 'up_in_pah']
    down_genes = degs[degs['regulation'] == 'down_in_pah']
    
    up_genes.to_csv(OUTPUT_DIR / "upregulated_in_pah.csv", index=False)
    down_genes.to_csv(OUTPUT_DIR / "downregulated_in_pah.csv", index=False)
    
    # Volcano plot
    plt.figure(figsize=(8, 6))
    colors = ['red' if r else 'gray' for r in mask]
    plt.scatter(result['log2fc_pah_vs_control'], -np.log10(0.05), c=colors, s=5, alpha=0.7)
    plt.xlabel('Log2 Fold Change (PAH vs Control)')
    plt.ylabel('-log10(p-value) [placeholder]')
    plt.title('Volcano Plot: PAH Lung DEGs (GSE117261)')
    plt.axvline(-1, color='k', linestyle='--', linewidth=0.5)
    plt.axvline(1, color='k', linestyle='--', linewidth=0.5)
    plt.savefig(OUTPUT_DIR / "pah_volcano_plot.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Step 2 complete. Results saved to {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()