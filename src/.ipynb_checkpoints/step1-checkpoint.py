#!/usr/bin/env python3
"""
Step 1: Identify lung- and heart-enriched genes from GTEx median TPM.

Input:
    data/processed/gtex_lung_heart_tpm.csv
        - gene_symbol
        - lung_tpm
        - heart_tpm

Output:
    results/step1_outputs/
        ├── differential_genes.csv
        ├── lung_enriched_genes.csv
        ├── heart_enriched_genes.csv
        ├── top20_lung_enriched.csv
        ├── top20_heart_enriched.csv
        └── lung_vs_heart_scatter.png
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Paths
INPUT_FILE = Path("data/processed/gtex_lung_heart_tpm.csv")
OUTPUT_DIR = Path("results/step1_outputs")
OUTPUT_DIR.mkdir(exist_ok=True)


def main():
    # Load processed GTEx data
    df = pd.read_csv(INPUT_FILE)
    
    # Validate columns
    required = ['gene_symbol', 'lung_tpm', 'heart_tpm']
    if not all(col in df.columns for col in required):
        raise ValueError(f"Input must contain columns: {required}")
    
    # Compute metrics
    df['mean_tpm'] = (df['lung_tpm'] + df['heart_tpm']) / 2.0
    df['log2fc'] = np.log2((df['lung_tpm'] + 1) / (df['heart_tpm'] + 1))
    
    # Filter differentially expressed genes
    mask = (df['mean_tpm'] > 1.0) & (np.abs(df['log2fc']) > 2.0)
    diff_genes = df[mask].copy()
    diff_genes['enrichment'] = diff_genes['log2fc'].apply(
        lambda x: 'lung_enriched' if x > 0 else 'heart_enriched'
    )
    
    # Save all DEGs
    diff_genes.to_csv(OUTPUT_DIR / "differential_genes.csv", index=False)
    
    # Split by tissue
    lung_genes = diff_genes[diff_genes['enrichment'] == 'lung_enriched']
    heart_genes = diff_genes[diff_genes['enrichment'] == 'heart_enriched']
    
    lung_genes.to_csv(OUTPUT_DIR / "lung_enriched_genes.csv", index=False)
    heart_genes.to_csv(OUTPUT_DIR / "heart_enriched_genes.csv", index=False)
    
    # Top 20 by |log2FC|
    lung_genes.nlargest(20, 'log2fc').to_csv(OUTPUT_DIR / "top20_lung_enriched.csv", index=False)
    heart_genes.nsmallest(20, 'log2fc').to_csv(OUTPUT_DIR / "top20_heart_enriched.csv", index=False)
    
    # Plot
    plt.figure(figsize=(6, 6))
    plt.scatter(df['lung_tpm'], df['heart_tpm'], alpha=0.5, s=1)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Lung TPM (log scale)')
    plt.ylabel('Heart TPM (log scale)')
    plt.title('Lung vs Heart Gene Expression (GTEx v8)')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.savefig(OUTPUT_DIR / "lung_vs_heart_scatter.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Step 1 complete. Results saved to {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()