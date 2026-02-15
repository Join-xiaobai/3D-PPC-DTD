#!/usr/bin/env python3
"""
Step 4: Compute lung vascular specificity score using GSE210248.

Input:
    data/processed/GSE210248_series_matrix.txt

Output:
    results/step4_outputs/
        └── lung_vascular_specificity_scores.csv
            - gene_symbol
            - vascular_mean (PAEC + PMVEC average)
            - whole_lung_mean
            - vascular_specificity_score (log2 ratio)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from io import StringIO

INPUT_FILE = Path("data/processed/GSE210248_series_matrix.txt")
OUTPUT_DIR = Path("results/step4_outputs")
OUTPUT_DIR.mkdir(exist_ok=True)


def parse_gse210248(filepath):
    """Parse GSE210248 to extract vascular vs whole lung expression."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Extract sample metadata
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
    
    # Load expression matrix
    data_start = next(i for i, line in enumerate(lines) if line.startswith('!series_matrix_table_begin')) + 1
    expr_lines = lines[data_start:-1]
    df = pd.read_csv(StringIO(''.join(expr_lines)), sep='\t', index_col=0)
    
    # Classify samples
    vascular_samples = []
    whole_lung_samples = []
    
    for col in df.columns:
        meta = samp_meta.get(col, {})
        tissue = meta.get('tissue', '').lower()
        cell_type = meta.get('cell type', '').lower()
        
        # Vascular: PAEC or PMVEC
        if 'pulmonary artery endothelial' in cell_type or 'lung microvascular endothelial' in cell_type:
            vascular_samples.append(col)
        # Whole lung
        elif 'lung' in tissue and 'whole' in tissue:
            whole_lung_samples.append(col)
    
    return df, vascular_samples, whole_lung_samples


def main():
    print("Loading GSE210248 (lung vascular vs whole lung)...")
    expr_df, vascular_samps, whole_samps = parse_gse210248(INPUT_FILE)
    
    if not vascular_samps or not whole_samps:
        raise ValueError("Could not identify vascular or whole lung samples in GSE210248.")
    
    vascular_mean = expr_df[vascular_samps].mean(axis=1)
    whole_lung_mean = expr_df[whole_samps].mean(axis=1)
    
    # Compute vascular specificity score: log2( (vascular + 1) / (whole_lung + 1) )
    specificity_score = np.log2((vascular_mean + 1) / (whole_lung_mean + 1))
    
    result = pd.DataFrame({
        'gene_symbol': vascular_mean.index,
        'vascular_mean': vascular_mean.values,
        'whole_lung_mean': whole_lung_mean.values,
        'vascular_specificity_score': specificity_score.values
    })
    
    # Save full score list
    result.to_csv(OUTPUT_DIR / "lung_vascular_specificity_scores.csv", index=False)
    
    # Top 100 vascular-enriched genes
    top_vascular = result.nlargest(100, 'vascular_specificity_score')
    top_vascular.to_csv(OUTPUT_DIR / "top100_lung_vascular_enriched.csv", index=False)
    
    # Plot distribution
    plt.figure(figsize=(8, 5))
    plt.hist(result['vascular_specificity_score'], bins=100, color='skyblue', edgecolor='k', alpha=0.7)
    plt.xlabel('Lung Vascular Specificity Score (log2 ratio)')
    plt.ylabel('Number of Genes')
    plt.title('Distribution of Lung Vascular Specificity Scores (GSE210248)')
    plt.axvline(1.0, color='red', linestyle='--', label='Threshold (2x enrichment)')
    plt.legend()
    plt.savefig(OUTPUT_DIR / "vascular_specificity_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Step 4 complete. Results saved to {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()