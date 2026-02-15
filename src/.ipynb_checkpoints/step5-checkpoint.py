#!/usr/bin/env python3
"""
Step 5: Integrate ChEMBL drug-target pairs with tissue- and disease-specific gene sets.

Inputs:
    data/processed/
        ├── chembl_drug_target_pairs.csv
        ├── uniprot_human_mapping.csv
    results/step1_outputs/lung_enriched_genes.csv
    results/step2_outputs/pah_lung_deg.csv
    results/step3_outputs/rv_cardiomyocyte_deg.csv
    results/step4_outputs/lung_vascular_specificity_scores.csv

Output:
    results/step5_outputs/
        └── candidate_drug_target_contexts.csv
            - drug_name, molecule_chembl_id
            - target_gene, target_chembl_id
            - pchembl_value
            - contexts: comma-separated list (e.g., "lung_enriched,pah_up,vascular_enriched")
"""

import pandas as pd
from pathlib import Path

# Paths
PROC = Path("data/processed")
STEP1 = Path("results/step1_outputs/lung_enriched_genes.csv")
STEP2 = Path("results/step2_outputs/pah_lung_deg.csv")
STEP3 = Path("results/step3_outputs/rv_cardiomyocyte_deg.csv")
STEP4 = Path("results/step4_outputs/lung_vascular_specificity_scores.csv")

OUTPUT_DIR = Path("results/step5_outputs")
OUTPUT_DIR.mkdir(exist_ok=True)


def main():
    print("Loading ChEMBL drug-target pairs...")
    chembl = pd.read_csv(PROC / "chembl_drug_target_pairs.csv")
    uniprot_map = pd.read_csv(PROC / "uniprot_human_mapping.csv")
    
    # Map target_chembl_id → UniProt → gene_symbol
    # Note: In ChEMBL v36, 'target_chembl_id' is already the target's ChEMBL ID,
    # but we need to link it to UniProt via target_dictionary (not done here for simplicity).
    # However, in many analyses, people use direct gene_symbol in activities.
    # Since our prep_data extracted only high-confidence pairs and assumed target_chembl_id maps to gene,
    # we will instead assume a simplified mapping: use a pre-linked file or note limitation.
    
    # ⚠️ WORKAROUND: In this pipeline, we assume that the target_chembl_id 
    # has already been mapped to gene_symbol during ChEMBL preprocessing.
    # But our current prep_data did NOT do that.
    # So we must revise: actually, we need to map via UniProt.
    
    # Let's correct: load full target dictionary to get component_uniprot
    # But to keep lightweight, we'll assume a common practice: 
    # Many public ChEMBL extracts include 'target_pref_name' as gene symbol.
    # Since we don't have that, we'll use an alternative: skip UniProt mapping for now,
    # and instead rely on the fact that in step5+ we only care about targets that appear in our gene lists.
    
    # Revised plan: We will not map ChEMBL targets to genes here.
    # Instead, we will intersect by gene symbol from our DEG lists.
    # But wait — our ChEMBL file has NO gene symbols!
    
    # 🔥 CRITICAL FIX: We must map ChEMBL targets to gene symbols.
    # Since we didn't do it in prep_data, we do minimal mapping here using UniProt.
    
    # Load full target dictionary (we have it in raw, but not processed)
    # To avoid complexity, let's assume we add a column in prep_data next time.
    # For now, we'll simulate by using a placeholder: we'll treat 'target_chembl_id' as proxy,
    # and only keep drugs whose targets are in our gene lists by name — which is impossible.
    
    # 💡 REALISTIC SOLUTION: In practice, ChEMBL activities often include 'target_pref_name'.
    # Let's adjust prep_data later. For now, we assume a corrected ChEMBL file that includes gene_symbol.
    
    # Given time, we'll assume that during ChEMBL processing, we already joined to gene_symbol.
    # So we modify expectation: `chembl_drug_target_pairs.csv` should contain 'target_gene_symbol'.
    
    # Since it doesn't, we abort and require fix.
    
    # BUT — to keep pipeline moving, we'll create a minimal mapping using UniProt:
    # Unfortunately, ChEMBL target_chembl_id → UniProt is in target_components table, which we didn't download.
    
    # 🛠 COMPROMISE: We will skip UniProt mapping and instead use the following:
    # Only consider targets that are **already in our gene expression matrices**,
    # and assume that the user has a way to link. For demo, we'll proceed by intersecting gene lists,
    # and assume ChEMBL targets are provided as gene symbols (common in simplified datasets).
    
    # Therefore, we assume `chembl_drug_target_pairs.csv` has a column 'target_gene_symbol'.
    # If not, this step fails.
    
    # Let's check columns
    if 'target_gene_symbol' not in chembl.columns:
        raise ValueError(
            "ChEMBL file must contain 'target_gene_symbol'. "
            "Please update prep_data.py to include UniProt mapping to gene symbols."
        )
    
    # Load gene sets
    lung_genes = set(pd.read_csv(STEP1)['gene_symbol'].astype(str))
    pah_deg = pd.read_csv(STEP2)
    pah_up = set(pah_deg[pah_deg['regulation'] == 'up_in_pah']['gene_symbol'].astype(str))
    pah_down = set(pah_deg[pah_deg['regulation'] == 'down_in_pah']['gene_symbol'].astype(str))
    
    rv_deg = pd.read_csv(STEP3)
    rv_up = set(rv_deg[rv_deg['regulation'] == 'up_in_pah_rv']['gene_symbol'].astype(str))
    rv_down = set(rv_deg[rv_deg['regulation'] == 'down_in_pah_rv']['gene_symbol'].astype(str))
    
    vascular = pd.read_csv(STEP4)
    vascular_enriched = set(vascular[vascular['vascular_specificity_score'] > 1.0]['gene_symbol'].astype(str))  # >2x
    
    # Annotate each drug-target pair with contexts
    def get_contexts(gene):
        ctx = []
        if gene in lung_genes:
            ctx.append('lung_enriched')
        if gene in pah_up:
            ctx.append('pah_lung_up')
        elif gene in pah_down:
            ctx.append('pah_lung_down')
        if gene in rv_up:
            ctx.append('pah_rv_up')
        elif gene in rv_down:
            ctx.append('pah_rv_down')
        if gene in vascular_enriched:
            ctx.append('lung_vascular_enriched')
        return ','.join(ctx) if ctx else 'other'
    
    chembl['contexts'] = chembl['target_gene_symbol'].astype(str).apply(get_contexts)
    candidates = chembl[chembl['contexts'] != 'other'].copy()
    
    # Save
    out_cols = [
        'drug_name', 'molecule_chembl_id',
        'target_gene_symbol', 'target_chembl_id',
        'pchembl_value', 'contexts'
    ]
    candidates[out_cols].to_csv(OUTPUT_DIR / "candidate_drug_target_contexts.csv", index=False)
    
    print(f"✅ Step 5 complete. Found {len(candidates)} context-relevant drug-target pairs.")
    print(f"Results saved to {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()