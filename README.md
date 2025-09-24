# PathFindeR Workflow

PathFindeR is an R-based workflow that expands a user supplied set of "seed" genes, constructs a protein interaction network, scores candidate genes with a random walk with restart (RWR), and performs multiple enrichment analyses to highlight putative disease modules. The repository bundles two scripts:

- **`PathFindeR_main.R`** — orchestrates the entire analysis pipeline.
- **`PathFindeR_Func.R`** — defines helper functions for working with STRING and NCBI resources, building graphs, running enrichment, and generating figures.

All results are written to `PF_output/`, including ranked gene tables, network exports, enrichment statistics, and figure files that can be explored or imported into external visualization tools.

## Required inputs

The main script expects several input files in the repository root. Replace these files with your own data (keeping the same structure) to reuse the workflow with a different disease context.

| File | Purpose | Key columns |
| ---- | ------- | ----------- |
| `My_genes_DecodeME_preprint.csv` | Curated list of 32 DecodeME candidate genes used as seeds. Semicolon separated. | `name`, `NCBI.id`, `list.name`, `weight`, `risk.locus`, `region` |
| `My_genes_DecodeME.csv` | Expanded fine-mapping catalogue from DecodeME/UK Biobank. Comma separated. Used to filter the preprint list and annotate supporting evidence. | `name`, `NCBI.id`, `list.name`, `weight`, `Gene_type`, `Study_type`, `Phenotype`, ... |
| `*-associated-targets-*.tsv` | Export from the Open Targets Platform for the disease of interest. Any file whose name contains `-associated-targets-` will be loaded. | `symbol` plus association metadata |
| `44321_2025_258_MOESM7_ESM.xlsx` *(optional)* | Proteomic results from external studies (UK Biobank in the bundled example). If present, sheets 2 and 3 are processed and saved as `Data/UKBB_proteomics.xlsx`. | `name`, statistical columns |

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
- Tissue enrichment summaries produced by `Tissue.ORA()`.

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
