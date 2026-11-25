# ADRB2–AF drug-target Mendelian randomization: analysis code and result tables

Replication materials for:

> Sun G, Shen Y. *The β2-adrenoceptor as a causal pathway to atrial fibrillation: a drug-target Mendelian randomization study informing long-acting β2-agonist prescribing in COPD.* In preparation for submission to the *European Journal of Clinical Pharmacology*.

## Contents

- `R_scripts/` — 14 R scripts: instrument extraction, harmonisation, primary MR, heterogeneity and pleiotropy, leave-one-out, single-SNP Wald ratios, the parallel ADRB1 analysis, forest-plot redrawing, and the numerical audit.
- `result_tables/` — 10 CSV files, split into `mr_af/` (ADRB2) and `mr_adrb1/` (ADRB1).
- `LICENSE`, `CITATION.cff`, `.zenodo.json`, `RELEASE_NOTES.md` — repository metadata.

## Target journal

European Journal of Clinical Pharmacology (Springer Nature; hybrid journal, subscription route, no article processing charge). The manuscript title is fixed as in `Short_Report_ADRB2_AF_MR_v7.docx`; the v8 manuscript will be reformatted to this journal's Instructions for Authors.

## Data sources

| Role | Dataset | Identifier | Details |
|---|---|---|---|
| Exposure | eQTLGen whole-blood cis-eQTL for ADRB2 | IEU OpenGWAS `eqtl-a-ENSG00000169252` | n = 31,684; 19,010 SNPs; GRCh37 |
| Outcome | FinnGen atrial fibrillation and flutter | IEU OpenGWAS `finn-b-I9_AF` | 22,068 cases / 116,926 controls |

Summary statistics were retrieved on 12 November 2026 through the IEU OpenGWAS API. Re-running the extraction scripts requires your own OpenGWAS JWT set as the `OPENGWAS_JWT` environment variable; **no token is included in this package**. All derived numbers are supplied as CSV, so the results can be checked without re-querying the API.

## Software

R 4.6.1; `TwoSampleMR` 0.7.13; `ieugwasr`; `ggplot2`.

The forest plots were drawn manually with `ggplot2`: `mr_forest_plot()` returns an empty figure when fewer than five instruments are available, which is the case here (three ADRB2 and four ADRB1 variants).

## Suggested reproduction order

1. `R_scripts/mr_af.R` → `mr_af2.R` → `mr_af3.R` → `mr_af4.R` — instrument selection, harmonisation, MR estimates.
2. `R_scripts/mr_sens.R` — heterogeneity, pleiotropy, leave-one-out, single-SNP Wald ratios.
3. `R_scripts/mr_plot.R` — scatter, forest and leave-one-out figures.
4. `R_scripts/adrb1_mr.R` … `adrb1_wald.R` — the parallel ADRB1 analysis.
5. `R_scripts/audit_mr.R` — recomputation of Cochran's Q and the MR-Egger intercept from the harmonised data.

## Citation and metadata

`CITATION.cff` enables GitHub's *Cite this repository* button. `.zenodo.json` controls the metadata of the Zenodo record generated from a GitHub Release; **Zenodo gives `.zenodo.json` precedence over `CITATION.cff`**, so keep the author list and title identical in both files. After the first Release, add the version DOI to both files and to the manuscript. The repository URL and the Zenodo version DOI together form the Data availability statement required in the journal's submission system (Springer Nature encourages a statement with a persistent link to archived data).

## Licence

Code: MIT. Derived result tables: CC0 1.0. The underlying summary statistics belong to eQTLGen and FinnGen and are governed by their own terms.
