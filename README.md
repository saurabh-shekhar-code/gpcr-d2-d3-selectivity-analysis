# Project 3: Deciphering GPCR Subtype Selectivity — Dopamine D2 vs D3 Binding Plasticity

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Target: GPCR D2/D3](https://img.shields.io/badge/Targets-D2_(6CM4)_%2F_D3_(3PBL)-purple.svg)](https://www.rcsb.org/structure/3PBL)
[![R^2: 0.963](https://img.shields.io/badge/Selectivity_R%C2%B2-0.963-brightgreen.svg)]()
[![Neuropharmacology](https://img.shields.io/badge/Domain-Neuropharmacology-blue.svg)]()

A structural pharmacology and computational CADD investigation into the **molecular determinants of subtype discrimination between the human Dopamine D2 and D3 receptors**.

---

## 🧠 Biological Context: The GPCR Selectivity Conundrum

Dopamine D2 and D3 receptors (Class A GPCRs) are high-value targets in schizophrenia, Parkinson's disease, and substance use disorders:
- **D2 Receptor Blockade:** Mediates classic antipsychotic efficacy in the dorsal striatum, but triggers severe **extrapyramidal motor side effects (EPS)**, tardive dyskinesia, and hyperprolactinemia.
- **D3 Receptor Selectivity:** Concentrated in the mesolimbic system (nucleus accumbens); selective D3 modulation improves cognitive and negative symptoms **without motor deficits**.
- **The Paradox:** The **orthosteric binding pocket (OBP)** is **100% conserved** across key dopamine-coordinating residues between D2 and D3 (Asp$^{3.32}$, Ser$^{5.42}$, Ser$^{5.46}$, Phe$^{6.51}$, His$^{6.55}$, Tyr$^{7.35}$).

Subtype selectivity can therefore **only** be achieved by exploiting structural plasticity in the **Secondary Binding Pocket (SBP)** and extracellular loop 2 (ECL2).

---

## 🔬 Key Computational Findings

### 1. Structural Hotspot Divergence
- **Hotspot 2.61:** **Val91 in D2** versus **Phe106 in D3**. The aromatic phenylalanine in D3 provides an extensive $\pi$-$\pi$ stacking and hydrophobic platform absent in D2.
- **ECL2 Lip:** Divergent gating residues (Ile183/Pro184 in D2 vs Leu181/Val180 in D3) alter pocket entrance volume.
- **Hotspot 7.39:** **Ile420 in D2** versus **Tyr377 in D3** creates an asymmetric polar hydrogen-bonding vector.

### 2. Energetics and Fingerprint Profiling (12 Ligands)
- **Non-Selective Antipsychotics (Haloperidol, Risperidone):** Bind primarily into the OBP; experience identical binding thermodynamics between D2 and D3 ($\Delta\Delta G \approx -0.5$ to $-0.9$ kcal/mol).
- **Bitopic Antagonists (SB-277011-A, S33084):** Span both the OBP and SBP, gaining massive stabilization from Phe106$^{2.61}$ and Tyr377$^{7.39}$ in D3 ($\Delta\Delta G \approx +2.9$ kcal/mol, $>100$-fold selectivity).
- **In Silico vs Experimental Correlation:** High correlation ($R^2 = 0.963$, Pearson $r = 0.981$) between computed $\Delta\Delta G$ and experimental $\log_{10}(K_{i, D2} / K_{i, D3})$.

---

## 📁 Repository Structure
```text
03_gpcr_selectivity_d2_d3/
├── GPCR_D2_vs_D3_Selectivity_Analysis.ipynb    # Interactive Jupyter Notebook
├── GPCR_D2_vs_D3_Selectivity_Analysis.html     # Standalone Interactive HTML Report
├── gpcr_selectivity_pipeline.py                # Reproducible analysis pipeline
├── build_project3_notebook_and_html.py          # Report and notebook builder
├── data/
│   ├── d2_d3_ligand_selectivity.csv            # 12 ligands x D2/D3 binding parameters
│   ├── d2_d3_pocket_alignment.csv              # Ballesteros-Weinstein sequence mapping
│   └── d2_d3_ifp_matrix.csv                    # Per-residue interaction fingerprints
└── figures/
    ├── d2_d3_selectivity_scatter.png           # Subtype selectivity scatter plot
    ├── d2_d3_sbp_interaction_fingerprint.png   # OBP vs SBP interaction heatmap
    ├── d2_d3_radar_profile.png                 # Multi-vector radar profile
    └── d2_d3_correlation_plot.png              # In silico vs in vitro correlation (R²=0.963)
```

---

## 🚀 Quick Reproduction

Run the complete pipeline:
```bash
python3 gpcr_selectivity_pipeline.py
```

Launch the interactive notebook:
```bash
jupyter notebook GPCR_D2_vs_D3_Selectivity_Analysis.ipynb
```
