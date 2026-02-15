#!/usr/bin/env python3
"""
Step 8: Functional enrichment analysis of top PAH candidate targets.

Input:
    results/step6_outputs/prioritized_pah_candidates.csv

Output:
    results/step8_outputs/
        ├── top100_pah_targets.txt                ← 基因列表
        ├── enrichment_results.csv                 ← 显著富集项（FDR < 0.05）
        └── enrichment_summary.txt                 ← 人类可读摘要
"""

import pandas as pd
from pathlib import Path
from gprofiler import GProfiler

INPUT_FILE = Path("results/step6_outputs/prioritized_pah_candidates.csv")
OUTPUT_DIR = Path("results/step8_outputs")
OUTPUT_DIR.mkdir(exist_ok=True)


def main():
    print("Loading top candidate targets...")
    df = pd.read_csv(INPUT_FILE)
    
    # Get top 100 unique target genes
    top_genes = df['target_gene_symbol'].drop_duplicates().head(100).tolist()
    
    # Save gene list
    gene_list_file = OUTPUT_DIR / "top100_pah_targets.txt"
    with open(gene_list_file, 'w') as f:
        f.write('\n'.join(top_genes))
    
    print(f"Running enrichment analysis on {len(top_genes)} genes...")
    
    # Initialize g:Profiler
    gp = GProfiler(return_dataframe=True)
    
    try:
        # Run enrichment: GO BP, KEGG, Reactome
        results = gp.profile(
            organism='hsapiens',
            query=top_genes,
            sources=['GO:BP', 'KEGG', 'REAC'],
            significance_threshold_method='fdr',
            no_iea=True,  # exclude electronic annotations
            max_p_value=1.0,
            min_set_size=5,
            max_set_size=500
        )
    except Exception as e:
        print(f"⚠️  Enrichment failed (no internet? g:Profiler down?): {e}")
        # Create empty placeholder
        results = pd.DataFrame(columns=[
            'source', 'term_id', 'term_name', 'p_value', 'significant'
        ])
    
    # Filter significant results (FDR < 0.05)
    if not results.empty:
        sig_results = results[results['significant']].copy()
        sig_results = sig_results.sort_values('p_value')
    else:
        sig_results = pd.DataFrame()
    
    # Save full enrichment results
    enrichment_file = OUTPUT_DIR / "enrichment_results.csv"
    sig_results.to_csv(enrichment_file, index=False)
    
    # Generate summary text
    summary_lines = []
    summary_lines.append("=" * 60)
    summary_lines.append("FUNCTIONAL ENRICHMENT OF TOP 100 PAH CANDIDATE TARGETS")
    summary_lines.append("=" * 60)
    summary_lines.append(f"Input genes: {len(top_genes)}")
    summary_lines.append(f"Significant terms (FDR < 0.05): {len(sig_results)}\n")
    
    if not sig_results.empty:
        # Group by source
        for source in ['GO:BP', 'KEGG', 'REAC']:
            sub = sig_results[sig_results['source'] == source].head(5)
            if not sub.empty:
                name_map = {'GO:BP': 'Gene Ontology (Biological Process)', 'KEGG': 'KEGG Pathways', 'REAC': 'Reactome Pathways'}
                summary_lines.append(f"Top enriched {name_map[source]}:")
                for _, row in sub.iterrows():
                    summary_lines.append(f"  • {row['term_name']} (FDR = {row['p_value']:.2e})")
                summary_lines.append("")
    else:
        summary_lines.append("No significant enrichment found (check gene list or connectivity).")
    
    # Write summary
    summary_file = OUTPUT_DIR / "enrichment_summary.txt"
    with open(summary_file, 'w') as f:
        f.write('\n'.join(summary_lines))
    
    print(f"✅ Step 8 complete.")
    print(f"   Gene list:      {gene_list_file.resolve()}")
    print(f"   Enrichment:     {enrichment_file.resolve()}")
    print(f"   Summary:        {summary_file.resolve()}")


if __name__ == "__main__":
    main()