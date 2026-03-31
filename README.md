# Spatial immune microenvironment in early-stage melanoma: D-ESMEL study

Analysis code accompanying the manuscript:

> **Spatial organization of the tumor immune microenvironment is associated 
> with distant metastatic outcome in stage I/II cutaneous melanoma**  
> Zhou et al. 

---

## Study overview

This is a matched case-control study (D-ESMEL, Dutch Early-Stage Melanoma Study) 
investigating the tumor immune microenvironment (TIME) of stage I/II primary 
cutaneous melanoma. Cases developed distant metastases; controls remained 
metastasis-free. Pairs were matched on Breslow thickness, ulceration, sex, and age.

Two analysis arms:
- **RNA-seq** (n = 178 complete pairs): bulk immune deconvolution (EPIC, CIBERSORT) 
  and TCR/BCR repertoire analysis (MiXCR)
- **Multiplex immunofluorescence (MxIF)** (n = 115 complete pairs): spatial analysis 
  of CD8, CD3, CD79a, CD68, and MelanA using SPIAT in R

---

## Repository structure

01_epic_cibersort.R         — Immune deconvolution analysis
02_tcr_bcr_repertoire.R     — TCR/BCR diversity analysis (MiXCR output)
03_spiat_density.R          — MxIF cell density analysis
04_spiat_mnnd.R             — Spatial distance (MNND) analysis
05_spiat_phenotypes.R       — Immune phenotype classification
06_hdbscan_clustering.R     — HDBSCAN cluster analysis

---

## Requirements

**R** (version 4.4.2)  
Key packages: `SPIAT`, `tidyverse`, `vegan`, `openxlsx`, `dbscan`

---

## Data availability

Patient-level data are not publicly available due to privacy regulations. 
Data may be available on reasonable request to the corresponding author 
(l.hollestein@erasmusmc.nl).

---

## Authors

Catherine Zhou, Thamila Kerkour, Alex Nigg, ,John Martens,
Karishma Lila, Hayri Balcioglu, Astrid Oostvogels, Reno Debets, Francesca Bosisio, 
Marcel Smid, Marlies Wakkee, Thierry P.P. van den Bosch, Antien Mooyaart, Loes Hollestein

Erasmus MC Cancer Institute, Rotterdam, the Netherlands
