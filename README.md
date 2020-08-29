# mtscATACpaper reproducibility repository

Updated: August 29, 2020

## About

This repository contains all code needed to reproduce all analyses [from this work](https://www.nature.com/articles/s41587-020-0645-6):


C. A. Lareau*, L. S. Ludwig*, _et al._ Massively parallel single-cell mitochondrial DNA genotyping and chromatin profiling. **Nature Biotechnology.** DOI: 10.1038/s41587-020-0645-6. 


Each folder corresponds to a particular dataset or analysis mode that is presented in this work. Rscripts can be found in the `code` folder that can be executed sequentially to reproduce analyses in part or in whole. Intermediate output files are also hosted for ease of use. 


## Setup

To best use this resource, we recommend pairing with large data files (that are not compatible with github as they exceed 100Mb). These files are available from the [Open Science Framework](https://osf.io/bupge/).

Once one downloads the zip archieve from OSF (~36 Gb), place the extracted folder in the same directory as this repository named `mtscATACpaper_large_data_files` (as shown below). This will enable running custom code to reproduce items in the output folders. 

```
.
├── mtscATACpaper_large_data_files
│   ├── intermediate
│   └── source
└── mtscATACpaper_reproducibility
    ├── cnv_compute
    ├── experimental_details
    ├── figure_CLL
    ├── figure_CRC
    ├── figure_invitro_CD34
    ├── figure_paired_cd34_pbmc
    ├── figures_TF1_GM11906_mixing
    ├── figure_supplement_public_scRNA
    ├── global_functions
    ├── numt_analysis
    └── README.md

```

The code assumes a relative file path with this organization. 


## See also

The main inputs are assumed to be output files from [mgatk](https://github.com/caleblareau/mgatk) and [CellRanger-ATAC](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac).

## Contact

For questions related to this work, please raise an issue on this repository. 

<br><br>
