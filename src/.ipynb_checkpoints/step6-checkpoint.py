#!/usr/bin/env python3
"""
Step 6: Score and prioritize drug-target candidates for PAH therapy.

Inputs:
    results/step5_outputs/candidate_drug_target_contexts.csv
    results/step1_outputs/lung_enriched_genes.csv
    results/step2_outputs/pah_lung_deg.csv
    results/step3_outputs/rv_cardiomyocyte_deg.csv
    results/step4_outputs/lung_vascular_specificity_scores.csv

Output:
    results/step6_outputs/
        └── prioritized_pah_candidates.csv
            - drug_name, molecule_chembl_id
            - target_gene_symbol, target_chembl_id
            - pchembl_value
            - lung_enriched (bool)
            - pah_lung_up (bool)
            - pah_rv_down (bool)      ← cardioprotective!
            - vascular_score (float)
            - composite_score (float)
"""

import pandas as pd
from pathlib import Path

STEP5 = Path("results/step5_outputs/candidate_drug_target_contexts.csv")
STEP1 = Path("results/step1_outputs/lung_enriched_genes.csv")
STEP2 = Path("results/step2_outputs/pah_lung_deg.csv")
STEP3 = Path("results/step3_outputs/rv_cardiomyocyte_deg.csv")
STEP4 = Path("results/step4_outputs/lung_vascular_specificity_scores.csv")

OUTPUT_DIR = Path("results/step6_outputs")
OUTPUT_DIR.mkdir(exist_ok=True)


def main():
    print("Loading candidate drug-target pairs...")
    candidates = pd.read_csv(STEP5)
    
    # Build gene annotation map
    gene_info = {}
    
    # Lung enriched (Step 1)
    lung_genes = set(pd.read_csv(STEP1)['gene_symbol'].astype(str))
    
    # PAH lung DEGs (Step 2)
    pah_lung = pd.read_csv(STEP2).set_index('gene_symbol')
    pah_lung_up = set(pah_lung[pah_lung['regulation'] == 'up_in_pah'].index.astype(str))
    
    # RV cardiomyocyte DEGs (Step 3)
    rv_deg = pd.read_csv(STEP3).set_index('gene_symbol')
    pah_rv_down = set(rv_deg[rv_deg['regulation'] == 'down_in_pah_rv'].index.astype(str))  # cardioprotective if upregulated
    
    # Vascular specificity (Step 4)
    vascular = pd.read_csv(STEP4).set_index('gene_symbol')['vascular_specificity_score'].to_dict()
    
    # Annotate each candidate
    def annotate_row(row):
        gene = str(row['target_gene_symbol'])
        return pd.Series({
            'lung_enriched': gene in lung_genes,
            'pah_lung_up': gene in pah_lung_up,
            'pah_rv_down': gene in pah_rv_down,
            'vascular_score': vascular.get(gene, 0.0)
        })
    
    annotations = candidates.apply(annotate_row, axis=1)
    df = pd.concat([candidates, annotations], axis=1)
    
    # Compute composite score
    # Rationale:
    # - Lung enrichment: +1
    # - PAH lung up (inhibitable): +2
    # - PAH RV down (activatable): +2 (cardioprotective!)
    # - Vascular score: add directly (log2 ratio)
    # - pChEMBL: already high (>6), used as weight later if needed
    df['composite_score'] = (
        df['lung_enriched'].astype(int) * 1 +
        df['pah_lung_up'].astype(int) * 2 +
        df['pah_rv_down'].astype(int) * 2 +
        df['vascular_score']
    )
    
    # Sort by composite score (descending)
    final = df.sort_values('composite_score', ascending=False)
    
    # Save full prioritized list
    out_cols = [
        'drug_name', 'molecule_chembl_id',
        'target_gene_symbol', 'target_chembl_id',
        'pchembl_value',
        'lung_enriched', 'pah_lung_up', 'pah_rv_down',
        'vascular_score', 'composite_score'
    ]
    final[out_cols].to_csv(OUTPUT_DIR / "prioritized_pah_candidates.csv", index=False)
    
    # Also save top 50
    final.head(50).to_csv(OUTPUT_DIR / "top50_pah_candidates.csv", index=False)
    
    print(f"✅ Step 6 complete. Top candidate:")
    if not final.empty:
        top = final.iloc[0]
        print(f"   Drug: {top['drug_name']} (ChEMBL: {top['molecule_chembl_id']})")
        print(f"   Target: {top['target_gene_symbol']} | Score: {top['composite_score']:.2f}")
    print(f"Results saved to {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()