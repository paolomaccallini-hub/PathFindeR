# PathFindeR Workflow

PathFindeR is an R-based workflow that expands a user supplied set of "seed" genes, constructs a protein interaction network, scores candidate genes with a random walk with restart (RWR), and performs multiple enrichment analyses to highlight putative disease modules. The repository bundles two scripts:

- **`PathFindeR_main.R`** — orchestrates the entire analysis pipeline.
- **`PathFindeR_Func.R`** — defines helper functions for working with STRING and NCBI resources, building graphs, running enrichment, and generating figures.

All results are written to `PF_output/`, including ranked gene tables, network exports, enrichment statistics, and figure files that can be explored or imported into external visualization tools.

## Required inputs

The main script expects several input files in the repository root. Replace these files with your own data (keeping the same structure) to reuse the workflow with a different disease context.

| File | Purpose | 
| ---- | ------- | 
| `My_genes_DecodeME_preprint.csv` | Curated list of 32 DecodeME candidate genes from the DecodeME preprint ([Ref](https://www.medrxiv.org/content/10.1101/2025.08.06.25333109v1)). | 
| `My_genes_DecodeME.csv` | Results of GWAS analysis and fine mapping of the DecodeME sumstats ([Ref](https://github.com/paolomaccallini-hub/DecodeME)). Used to filter the preprint list and annotate supporting evidence. |
| `OT-EFO_0004540-associated-targets-7_30_2025-v25_06.tsv` | Export from the Open Targets Platform for the ME/CFS. Any file whose name contains `-associated-targets-` will be loaded. |
| `44321_2025_258_MOESM7_ESM.xlsx` | Proteomic results from UK Biobank ([Ref](https://www.embopress.org/doi/full/10.1038/s44321-025-00258-8)). If present, sheets 2 and 3 are processed and saved as `Data/UKBB_proteomics.xlsx`. | 
| `media-9.xlsx` | A disease module for ME/CFS of 115 genes, from a recent publication on 400 patients ([Ref](https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/)).  |

### External databases

During its first run the pipeline downloads auxiliary datasets automatically and caches them under `Data/`:

- STRING v12 protein links (`Data/STRING/9606.protein.links.v12.0.txt.gz`) and protein annotations.
- NCBI human gene information (`Data/NCBI/human_gene_info.gz`).

No manual preparation is needed, but the first execution will take longer while these resources are fetched.

## Software prerequisites

Install the required R packages before running the workflow. The scripts rely on commonly available Bioconductor and CRAN libraries, including:

`dplyr`, `rentrez`, `httr`, `jsonlite`, `curl`, `biomaRt`, `stringr`, `DOSE`, `ReactomePA`, `igraph`, `MASS`, `data.table`, `clusterProfiler`, `org.Hs.eg.db`, `pathview`, `enrichplot`, `GOplot`, `readxl`, `writexl`, and `TissueEnrich`.

Bioconductor packages can be installed with `BiocManager::install(<package_name>)`, while CRAN packages can be installed with `install.packages()`.

The scripts also assume access to the internet to download STRING and NCBI databases and, optionally, to resolve gene identifiers via the NCBI Entrez API.

## Running the pipeline

1. Clone or copy this repository and place the required input files in the project root (see table above).
2. Ensure all packages listed in the prerequisites section are installed.
3. From the repository root, run the main script:

   ```bash
   Rscript PathFindeR_main.R
   ```

The script sources `PathFindeR_Func.R`, creates any missing subdirectories under `PF_output/`, builds the STRING interaction network, aligns input identifiers to STRING preferred names, and expands the seed set.

## What the script produces

Running `PathFindeR_main.R` generates a comprehensive set of outputs inside `PF_output/`:

- `All_genes.tsv` — every gene discovered during expansion with RWR scores, ranks, cumulative stationary probabilities, and provenance tags indicating whether a gene originated as a seed or was predicted via STRING.
- `Disease_results_list.rds` — serialized list containing the curated seed table, expanded gene table, adjacency matrix, and expanded gene-list matrix for downstream reuse.
- Graph exports (`All_genes_graph.tiff`, `All_genes_cytoscape.tsv`, `Disease_module_graph.tiff`, etc.) describing the merged network and the high-confidence disease module.
- ORA/GSEA outputs under `PF_output/ORA/` and `PF_output/GSEA/`, including tables of enriched KEGG, Reactome, GO, and Disease Ontology pathways and corresponding visualization files. KEGG pathway diagrams are moved into `PF_output/ORA/KEGG/`, while gene lists ready for Reactome’s Pathway Browser are saved into `PF_output/ORA/Reactome/`.

## Customising the analysis

- **Different seeds**: replace the contents of `My_genes_DecodeME_preprint.csv` and (optionally) `My_genes_DecodeME.csv` with your own gene lists. Maintain the column headers so the script can read them without modification.
- **Alternative evidence sources**: swap in another Open Targets export or omit it entirely (the script simply skips associations not present in STRING).
- **Scoring thresholds**: parameters such as the cumulative stationary probability cut-offs (`0.80`, `0.85`, `0.90`, `0.95`) are defined in `PathFindeR_main.R` and can be edited to change how stringent the disease module selection is.
- **Enrichment depth**: adjust `top.results` and `list.number` to report more (or fewer) pathways per database for GSEA and ORA.

The random seed is fixed (`set.seed(12345)`) to guarantee reproducible results across runs with the same inputs and software versions.

## Troubleshooting tips

- If identifier lookups fail with HTTP 500 errors, the helper functions automatically retry up to 100 times. Persistent failures usually indicate a temporary NCBI outage.
- When introducing new genes, confirm that they have valid NCBI IDs and STRING preferred names—messages in the R console will highlight entries that cannot be resolved.
- Large STRING downloads can take several minutes. Check your network connection if the script appears idle during the initial setup.

For any changes to the workflow, edit `PathFindeR_main.R` (for orchestration logic) or `PathFindeR_Func.R` (for helper routines) and rerun the script to regenerate outputs.

## Methods and Results

### Gene expansion, adjacency matrix, and graph

All the genes with a PPI score above or equal to 0.7 are retrieved from a local installation of the STRING database (v12) for the 18 seeds (RABGAP1L, DARS2, RC3H1, ZBTB37, TNFSF4, ANKRD45, KLHL20, PRDX6, SERPINC1, SLC9C2, BTN2A2, TRIM38, ABT1, OLFM4, CSE1L, ARFGEF2, STAU1, ZNFX1). This leads to a merged gene list (MGL) of 345 genes, available in this repo as `All_genes.tsv`. For each predicted gene, this table contains the highest PPI score. It also keeps track of the seeds that generated each predicted gene, and for each seed, the complete set of predicted genes is indicated. After gene expansion, all the missing PPI scores are retrieved and stored in a symmetric matrix (adjacency matrix), saved as `Gene_matrix.tsv`, also present in this repository. A plot of the graph associated with the adjacency matrix is generated and reported below. Graphs are generated using the package `igraph`, and a file to input into Cytoscape is also generated (available in this repository as `All_genes_cytoscape.tsv`). 

![All_genes_graph](https://github.com/user-attachments/assets/f521eca0-5034-4097-b478-f2f66ab9cf8a)

In the image above, note that: seeds are written with a red font; genes that are specific to the GWAS obtained using patients who reported an infection at the onset are indicated in green, while genes associated with GWAS-1 (the main cohort of DecodeME) are indicated in purple.

### Random walk with restart and pruning

Genes of the MGL are ranked by the stationary probability of a random walk with restart, with initial probability assigned only to the seeds. Note that seeds that were ranked as tier one in the DecodeME preprint were assigned an initial probability double that of the one assigned to genes ranked tier two. In our case, only one of the 18 genes selected by fine-mapping is a tier two gene (OLFM4). The back probability employed is 0.7, and the stationary probability is calculated by the function `page_rank`. Then, genes are ranked according to their stationary probability (see column `score.RWR` of file `All_genes.tsv`). By selecting 90% of the cumulative stationary probability (CSP), a disease module of 119 high-ranking genes is selected and plotted (see figure below). The corresponding files (as for the full MGL) are available in this repo.

![Disease_module_graph](https://github.com/user-attachments/assets/b3e6f98c-7f83-407e-bccc-7e32613702e9)

### Gene-set enrichment analysis

The script performs gene-set enrichment analysis (GSEA) on the MGL ranked by RWR. It uses the STRING database as a background and works on Reactome, KEGG, Disease Ontology, and Gene Ontology. No significant enrichment is reported.

### Over-representation analysis

Over-representation analysis on Reactome, KEGG, Disease Ontology, and Gene Ontology is performed on the disease module. The top results are available in this repo, file `ORA_RWR.tsv`. The script also generates several plots to summarise the main results, which are reported below.

![KEGG enrichment](https://github.com/user-attachments/assets/ff09f61a-098f-49dd-8dcf-43957a65d425)
![Reactome enrichment](https://github.com/user-attachments/assets/35d78f83-312a-4516-952b-05ccfdf22960)
![DO enrichment](https://github.com/user-attachments/assets/2a2ed630-081f-4692-a2ea-3833ab537cdc)
![GO enrichment](https://github.com/user-attachments/assets/cc3174ab-c285-4a7d-b5f4-76afb234d351)
![GO_ORA_goplot](https://github.com/user-attachments/assets/8c218bd4-28e9-475f-989b-8a236612f30e)
![GO_ORA_cnetplot](https://github.com/user-attachments/assets/d59d99e0-5198-40b3-8b30-a6b4a81ec7e6)

### KEGG pathways

The top-ten KEGG pathways with a significant overlap from over-representation analysis are plotted with the enzymes from the disease module highlighted in red. Plotting are available in the file section of this repo. I've included one of them below.

<img width="1072" height="895" alt="hsa03013 Disease_ORA" src="https://github.com/user-attachments/assets/b5038f0f-36b1-4309-8b2d-87ff7a2ca573" />

### Comparison with results from previous experiments

I compared the MGL and disease modules selected at different cut-offs of CSP (80%, 85%, 90%, 95%) with the following results from previous experiments:

| Description                                       | Input File | Reference |
| ------------------------------------------------- | ---------- | --------- | 
| Disease module of 115 genes built using whole-genome sequencing data from 400 ME/CFS patients | `media-9.xlsx` | [Ref](https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/) |
| Proteomic study on 1455 ME/CFS cases from the UK Biobank database, compared with 131,303 controls. I used the total sample (males plus females) and the total effect |  `44321_2025_258_MOESM7_ESM.xlsx` | [Ref](https://www.embopress.org/doi/full/10.1038/s44321-025-00258-8) |
| Collection of 497 known disease-gene associations for ME/CFS from different experiments, according to Open Target Platform | `OT-EFO_0004540-associated-targets-7_30_2025-v25_06.tsv` | [Ref](https://platform.opentargets.org/)

Comparisons were performed using the hypergeometric test, using as background a set of 15,000 genes for the first and the third case, and of 2,895 (number of proteins tested in the proteomic study) for the second one. The results are collected in the file `replication.tsv` available in this repo and reported below:

| setA               | setB        | sizeA | sizeB | overlap | background | pvalue | genes                                    |
| ------------------ | ----------- | ----- | ----- | ------- | ---------- | ------ | ---------------------------------------- |
| All genes | Zhang S 2025 | 345 | 115 | 1 | 15000 | 9.32e-01 | CDC23 |
| Disease.module.0.8 | Zhang S 2025 | 23 | 115 | 0 | 15000 | 1.00e+00 | |
| Disease.module.0.85 | Zhang S 2025 | 50 | 115 | 0 | 15000 | 1.00e+00 | |
| Disease.module.0.9 | Zhang S 2025 | 119 | 115 | 0 | 15000 | 1.00e+00 | |
| Disease.module.0.95 | Zhang S 2025 | 212 | 115 | 0 | 15000 | 1.00e+00 | |
| All genes | Open Targets | 345 | 497 | 21 | 15000 | 5.59e-03 | ALB/CD28/CD4/CD40LG/CRP/CTLA4/ELANE/FLNA/FMR1/GSR/HPGDS/HPX/HSPA8/IL10/IL6/LTF/MTHFR/SERPINA5/SOD2/TNF/VARS2 |
| Disease.module.0.8 | Open Targets | 23 | 497 | 0 | 15000 | 1.00e+00 | |
| Disease.module.0.85 | Open Targets | 50 | 497 | 4 | 15000 | 8.32e-02 | FLNA/HSPA8/LTF/TNF |
| Disease.module.0.9 | Open Targets | 119 | 497 | 12 | 15000 | 5.84e-04 | CD28/CD4/CD40LG/CTLA4/FLNA/FMR1/GSR/HSPA8/IL10/IL6/LTF/TNF |
| Disease.module.0.95 | Open Targets | 212 | 497 | 14 | 15000 | 1.12e-02 | ALB/CD28/CD4/CD40LG/CTLA4/FLNA/FMR1/GSR/HSPA8/IL10/IL6/LTF/TNF/VARS2 |
| All genes | Proteomics | 109 | 233 | 17 | 2895 | 5.40e-03 | AGXT/AMBP/CD4/CD80/CHGA/F7/FGA/GGT1/GSTA1/GSTA3/HRG/PLAT/PROC/PROCR/SERPIND1/SERPINF2/TNFRSF9 |
| Disease.module.0.8 | Proteomics | 8 | 233 | 0 | 2895 | 1.00e+00 | |
| Disease.module.0.85 | Proteomics | 12 | 233 | 0 | 2895 | 1.00e+00 | |
| Disease.module.0.9 | Proteomics | 34 | 233 | 3 | 2895 | 5.24e-01 | CD4/CHGA/TNFRSF9 |
| Disease.module.0.95 | Proteomics | 55 | 233 | 5 | 2895 | 4.58e-01 | CD4/CD80/CHGA/FGA/TNFRSF9 |
