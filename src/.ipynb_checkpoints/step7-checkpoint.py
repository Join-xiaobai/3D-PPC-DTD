#!/usr/bin/env python3
"""
Step 7: Generate final candidate table and a plain-text summary report.

Input:
    results/step6_outputs/prioritized_pah_candidates.csv

Outputs:
    results/step7_outputs/
        ├── final_pah_candidate_list.csv          ← 完整候选表（含注释）
        └── pah_target_summary.txt                ← 人类可读的摘要
"""

import pandas as pd
from pathlib import Path

INPUT_FILE = Path("results/step6_outputs/prioritized_pah_candidates.csv")
OUTPUT_DIR = Path("results/step7_outputs")
OUTPUT_DIR.mkdir(exist_ok=True)


def main():
    print("Generating final candidate list and summary...")
    df = pd.read_csv(INPUT_FILE)

    # Add therapeutic rationale column
    def rationale(row):
        parts = []
        if row['pah_lung_up']:
            parts.append("Inhibit (↑ in PAH lung → pathogenic)")
        if row['pah_rv_down']:
            parts.append("Activate (↓ in PAH RV → cardioprotective loss)")
        if not parts:
            parts.append("Context unclear")
        return "; ".join(parts)

    df['therapeutic_rationale'] = df.apply(rationale, axis=1)

    # Reorder for clarity
    out_cols = [
        'drug_name', 'molecule_chembl_id',
        'target_gene_symbol', 'target_chembl_id',
        'pchembl_value',
        'lung_enriched', 'pah_lung_up', 'pah_rv_down',
        'vascular_score', 'composite_score',
        'therapeutic_rationale'
    ]
    final_df = df[out_cols].copy()

    # Save full candidate list
    output_csv = OUTPUT_DIR / "final_pah_candidate_list.csv"
    final_df.to_csv(output_csv, index=False)

    # Generate plain-text summary
    top5 = final_df.head(5)
    summary = []
    summary.append("=" * 60)
    summary.append("PULMONARY ARTERIAL HYPERTENSION (PAH) TARGET PRIORITIZATION")
    summary.append("=" * 60)
    summary.append(f"Total candidates with known drugs: {len(final_df)}")
    summary.append(f"Top 5 high-priority targets:\n")

    for i, (_, row) in enumerate(top5.iterrows(), 1):
        summary.append(f"{i}. Drug: {row['drug_name']} (ChEMBL: {row['molecule_chembl_id']})")
        summary.append(f"   Target: {row['target_gene_symbol']}")
        summary.append(f"   Rationale: {row['therapeutic_rationale']}")
        summary.append(f"   Composite Score: {row['composite_score']:.2f} | pChEMBL: {row['pchembl_value']}")
        summary.append("")

    summary.append("Scoring components:")
    summary.append("- Lung enriched: expressed higher in lung vs heart (GTEx)")
    summary.append("- PAH lung up: upregulated in PAH lung tissue (GSE117261)")
    summary.append("- PAH RV down: downregulated in PAH right ventricle cardiomyocytes (GSE240921)")
    summary.append("- Vascular score: log2(PAEC+PMVEC / whole lung) from GSE210248")
    summary.append("\nInterpretation:")
    summary.append("• 'Inhibit' targets are potential drivers of PAH pathology.")
    summary.append("• 'Activate' targets may restore cardioprotective function in the right heart.")

    # Write summary
    summary_file = OUTPUT_DIR / "pah_target_summary.txt"
    with open(summary_file, 'w') as f:
        f.write('\n'.join(summary))

    print(f"✅ Step 7 complete.")
    print(f"   Final candidate list: {output_csv.resolve()}")
    print(f"   Summary report:       {summary_file.resolve()}")


if __name__ == "__main__":
    main()