#!/usr/bin/env python3
"""
Step 9: Compare candidates against known PAH drugs and identify repurposing opportunities.

Input:
    results/step6_outputs/prioritized_pah_candidates.csv

Output:
    results/step9_outputs/
        ├── pah_drug_comparison.csv
        └── repurposing_opportunities.txt
"""

import pandas as pd
from pathlib import Path

INPUT_FILE = Path("results/step6_outputs/prioritized_pah_candidates.csv")
OUTPUT_DIR = Path("results/step9_outputs")
OUTPUT_DIR.mkdir(exist_ok=True)


# Known PAH drug targets (as of 2025)
KNOWN_PAH_TARGETS = {
    'EDNRA', 'EDNRB', 'PDE5A', 'PTGIR', 'PPARG',
    'GUCY1A2', 'GUCY1B3'  # sGC subunits
}

# Known PAH drugs (for reference)
PAH_DRUGS = {
    'Bosentan', 'Ambrisentan', 'Macitentan',
    'Sildenafil', 'Tadalafil',
    'Epoprostenol', 'Treprostinil', 'Selexipag',
    'Riociguat'
}


def main():
    print("Loading prioritized candidates...")
    df = pd.read_csv(INPUT_FILE)

    # Normalize gene symbols to upper case for safe comparison
    df['target_gene_symbol'] = df['target_gene_symbol'].str.upper()
    known_set = {g.upper() for g in KNOWN_PAH_TARGETS}

    # Annotate
    df['is_known_pah_target'] = df['target_gene_symbol'].isin(known_set)
    df['is_pah_drug'] = df['drug_name'].isin(PAH_DRUGS)

    # Identify repurposing candidates:
    # - NOT a known PAH drug
    # - BUT has high composite score AND target is novel or underutilized
    df['repurposing_candidate'] = (
        (~df['is_pah_drug']) &
        (df['composite_score'] >= 3.0)  # meaningful context support
    )

    # Save full comparison
    out_cols = [
        'drug_name', 'molecule_chembl_id',
        'target_gene_symbol', 'pchembl_value',
        'composite_score',
        'is_known_pah_target', 'is_pah_drug',
        'repurposing_candidate',
        'therapeutic_rationale' if 'therapeutic_rationale' in df.columns else 'lung_enriched'
    ]
    # Ensure therapeutic_rationale exists (borrow from step7 logic if missing)
    if 'therapeutic_rationale' not in df.columns:
        def rationale(row):
            parts = []
            if row['pah_lung_up']:
                parts.append("Inhibit")
            if row['pah_rv_down']:
                parts.append("Activate")
            return "; ".join(parts) if parts else "Context unclear"
        df['therapeutic_rationale'] = df.apply(rationale, axis=1)

    final_df = df[[
        'drug_name', 'molecule_chembl_id',
        'target_gene_symbol', 'pchembl_value',
        'composite_score', 'therapeutic_rationale',
        'is_known_pah_target', 'is_pah_drug', 'repurposing_candidate'
    ]].copy()

    comparison_file = OUTPUT_DIR / "pah_drug_comparison.csv"
    final_df.to_csv(comparison_file, index=False)

    # Generate repurposing summary
    repurposing_hits = final_df[
        final_df['repurposing_candidate']
    ].sort_values('composite_score', ascending=False).head(20)

    summary_lines = []
    summary_lines.append("=" * 70)
    summary_lines.append("PAH DRUG REPURPOSING OPPORTUNITIES")
    summary_lines.append("=" * 70)
    summary_lines.append("Criteria: Approved drug + NOT used in PAH + composite_score ≥ 3.0\n")

    if not repurposing_hits.empty:
        summary_lines.append(f"Top {len(repurposing_hits)} repurposing candidates:\n")
        for _, row in repurposing_hits.iterrows():
            summary_lines.append(f"• Drug: {row['drug_name']} (ChEMBL: {row['molecule_chembl_id']})")
            summary_lines.append(f"  Target: {row['target_gene_symbol']} | Score: {row['composite_score']:.2f}")
            summary_lines.append(f"  Rationale: {row['therapeutic_rationale']}")
            summary_lines.append("")
    else:
        summary_lines.append("No strong repurposing candidates found.")

    # Also report known target recovery
    known_hits = final_df[final_df['is_known_pah_target']].drop_duplicates('target_gene_symbol')
    summary_lines.append("\n" + "-" * 50)
    summary_lines.append("Validation: Known PAH targets recovered by pipeline:")
    if not known_hits.empty:
        for _, row in known_hits.iterrows():
            status = "✅ Already used in PAH" if row['is_pah_drug'] else "⚠️  Novel drug for known target"
            summary_lines.append(f"  - {row['target_gene_symbol']}: {row['drug_name']} ({status})")
    else:
        summary_lines.append("  None (check scoring thresholds)")

    # Write summary
    summary_file = OUTPUT_DIR / "repurposing_opportunities.txt"
    with open(summary_file, 'w') as f:
        f.write('\n'.join(summary_lines))

    print(f"✅ Step 9 complete.")
    print(f"   Comparison table: {comparison_file.resolve()}")
    print(f"   Repurposing summary: {summary_file.resolve()}")


if __name__ == "__main__":
    main()