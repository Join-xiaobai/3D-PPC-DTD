# A Tissue-Selective Dual-Axis Therapeutic Strategy for Pulmonary Arterial Hypertension Enabled by 3D Pharmacological Point Cloud Modeling

This repository implements a **computational pipeline** to identify high-confidence, tissue-contextualized drug targets for **Pulmonary Arterial Hypertension (PAH)** by integrating multi-omics data across three critical biological axes:

1. **Lung vs. Heart Specificity** (GTEx)  
2. **Disease Dysregulation in PAH Lung** (GSE117261)  
3. **Cardioprotective Signaling in Right Ventricle** (GSE240921 + GSE210248)

We model each drug–target pair as a point in a **3D pharmacological point cloud**, where coordinates reflect:
- **X**: Pulmonary enrichment  
- **Y**: Pathogenic upregulation in PAH lung  
- **Z**: Cardioprotective downregulation in PAH right ventricle  

The result is a **dual-axis therapeutic strategy**:  
✅ **Inhibit** lung-enriched, PAH-upregulated drivers  
❤️ **Activate** cardioprotective pathways lost in right heart failure  

We further prioritize candidates with **lung vascular specificity** and **existing clinical compounds** (ChEMBL), enabling rapid repurposing.

---

## 📁 Project Structure

```
3D-PPC-DTD/
├── ipynb/ * The notbook version of the project is for reference only. Please see python code for details.
├── data/
│   ├── raw/
│   │   ├── GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm.gct.gz
│   │   ├── GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt
│   │   ├── chembl_36_activities.csv.gz
│   │   ├── chembl_36_molecule_dictionary.csv.gz
│   │   ├── chembl_36_target_dictionary.csv.gz
│   │   ├── chembl_36_target_components.csv.gz
│   │   ├── uniprot_human.tsv.gz
│   │   ├── GSE117261_series_matrix.txt.gz
│   │   ├── GSE240921_series_matrix.txt.gz
│   │   └── GSE210248_series_matrix.txt.gz
│   │
│   └── processed/
│       ├── gtex_lung_heart_tpm.csv
│       ├── GSE117261_series_matrix.txt
│       ├── GSE240921_series_matrix.txt
│       ├── GSE210248_series_matrix.txt
│       ├── chembl_drug_target_pairs.csv
│       └── uniprot_human_mapping.csv
│
├── src/
│   ├── prep_data.py
│   ├── step1.py
│   ├── step2.py
│   ├── step3.py
│   ├── step4.py
│   ├── step5.py
│   ├── step6.py
│   ├── step7.py
│   ├── step8.py
│   └── step9.py
│
└── results/
    ├── step1_outputs/
    │   ├── differential_genes.csv
    │   ├── lung_enriched_genes.csv
    │   ├── heart_enriched_genes.csv
    │   ├── top20_lung_enriched.csv
    │   ├── top20_heart_enriched.csv
    │   └── lung_vs_heart_scatter.png
    │
    ├── step2_outputs/
    │   ├── pah_lung_deg.csv
    │   ├── upregulated_in_pah.csv
    │   ├── downregulated_in_pah.csv
    │   └── pah_volcano_plot.png
    │
    ├── step3_outputs/
    │   ├── rv_cardiomyocyte_deg.csv
    │   ├── up_in_pah_rv.csv
    │   ├── down_in_pah_rv.csv
    │   └── rv_cardiomyocyte_deg_plot.png
    │
    ├── step4_outputs/
    │   ├── lung_vascular_specificity_scores.csv
    │   ├── top100_lung_vascular_enriched.csv
    │   └── vascular_specificity_distribution.png
    │
    ├── step5_outputs/
    │   └── candidate_drug_target_contexts.csv
    │
    ├── step6_outputs/
    │   ├── prioritized_pah_candidates.csv
    │   └── top50_pah_candidates.csv
    │
    ├── step7_outputs/
    │   ├── final_pah_candidate_list.csv
    │   └── pah_target_summary.txt
    │
    ├── step8_outputs/
    │   ├── top100_pah_targets.txt
    │   ├── enrichment_results.csv
    │   └── enrichment_summary.txt
    │
    └── step9_outputs/
        ├── pah_drug_comparison.csv
        └── repurposing_opportunities.txt
```


> ✅ All outputs are **plain-text or CSV** — no HTML, no databases. Fully reproducible and publication-ready.

---

## 🧪 Pipeline Overview

| Step | Function | Input(s) | Output(s) |
|------|--------|--------|---------|
| **Step 1** | Lung/Heart enrichment (GTEx) | `gtex_lung_heart_tpm.csv` | `lung_enriched_genes.csv`, scatter plot |
| **Step 2** | PAH lung DEGs | `GSE117261_series_matrix.txt` | `pah_lung_deg.csv`, volcano plot |
| **Step 3** | RV cardiomyocyte DEGs | `GSE240921_series_matrix.txt` | `rv_cardiomyocyte_deg.csv` |
| **Step 4** | Lung vascular specificity | `GSE210248_series_matrix.txt` | `lung_vascular_specificity_scores.csv` |
| **Step 5** | Map ChEMBL drugs to genes | `chembl_drug_target_pairs.csv` | `candidate_drug_target_contexts.csv` |
| **Step 6** | **3D Point Cloud Scoring** | Steps 1–5 outputs | `prioritized_pah_candidates.csv` |
| **Step 7** | Final candidate table + summary | Step 6 output | `final_pah_candidate_list.csv`, `pah_target_summary.txt` |
| **Step 8** | Functional enrichment (g:Profiler) | Top 100 targets | `enrichment_results.csv`, `enrichment_summary.txt` |
| **Step 9** | Clinical validation & repurposing | Step 6 output | `pah_drug_comparison.csv`, `repurposing_opportunities.txt` |

> 🔑 **Core innovation**: The **composite score** in Step 6 integrates all axes into a single prioritization metric:
> ```
> Score = LungEnriched×1 + PAHLungUp×2 + PAHRVDown×2 + VascularScore
> ```

---

## 🚀 Quick Start

### Prerequisites
- Python 3.8 or higher
- Internet access (required only for Step 8: g:Profiler API)

### Installation

```bash
git clone https://github.com/yourname/pah-point-cloud.git
cd pah-point-cloud
pip install -r requirements.txt
```

## Prepare Data

Manually download the following public datasets and place them in `data/raw/`:

| Dataset                  | Source                                                                 |
|--------------------------|------------------------------------------------------------------------|
| GTEx v8 TPM + annotations| [GTEx Portal](https://gtexportal.org/)                                |
| GSE117261                | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117261)   |
| GSE240921                | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240921)   |
| GSE210248                | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210248)   |
| ChEMBL v36               | [ChEMBL FTP](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_36/) |
| UniProt Human            | [UniProt](https://www.uniprot.org/proteomes/UP000005640)              |

Then run preprocessing:
```bash
python src/prep_data.py
```

## Run Full Analysis

```bash
python src/step1.py
python src/step2.py
python src/step3.py
python src/step4.py
python src/step5.py
python src/step6.py
python src/step7.py
python src/step8.py   # requires internet
python src/step9.py
```

⏱️ Total runtime: ~15–30 minutes on a standard laptop.

## 📊 Key Outputs to Explore

After completion, inspect these files:

- **`results/step6_outputs/prioritized_pah_candidates.csv`**  
  → Ranked list of drug–target pairs with composite scores.

- **`results/step7_outputs/pah_target_summary.txt`**  
  → Human-readable top-5 summary with therapeutic rationale.

- **`results/step9_outputs/repurposing_opportunities.txt`**  
  → High-potential "old drug, new use" candidates (e.g., kinase inhibitors, antifibrotics).

- **`results/step8_outputs/enrichment_summary.txt`**  
  → Reveals enriched pathways (e.g., *vascular smooth muscle contraction*, *response to hypoxia*).

## 📦 Dependencies (`requirements.txt`)

```txt
pandas>=1.5.0
numpy>=1.21.0
scipy>=1.9.0
matplotlib>=3.6.0
gprofiler-official>=1.0.0
```

```bash
pip install -r requirements.txt
```
