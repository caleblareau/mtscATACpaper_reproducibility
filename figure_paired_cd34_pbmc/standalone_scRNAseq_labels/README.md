
### Step 1. Raw data obtained from 10x website

```
https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3
```

### Step 2. Scrublet was used to remove multiplets

```
python scrublet_10x.py pbmc_10k_v3
```

### Step 3. Seurat / Rscript used to find and annotate clusters; produce RDS for scATAC integration

```
Rscript 01_seurat_analysis.R
```