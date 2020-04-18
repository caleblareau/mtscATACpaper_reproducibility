# mtscATACpaper reproducibility repository

Updated: April 18, 2020

## About

This repository contains all code needed to reproduce all analyses from this work:


C. A. Lareau*, L. S. Ludwig*, _et al._ Massively parallel joint single-cell mitochondrial DNA genotyping and chromatin profiling reveals properties of human clonal variation. **Under review.** 2020. 


Each folder corresponds to a particular dataset or analysis mode that is presented in this work. Rscripts can be found in the `code` folder that can be executed sequentially to reproduce analyses in part or in whole. Intermediate output files are also hosted for ease of use. 


## Setup

To best use this resource, we recommend pairing with large data files (that are not compatible with github as they exceed 100Mb). A `.tar.gz` of these files is [available at this link](<<dropbox link>>).

Once this file is downloaded (~36 Gb), place the extracted folder in the same directory as this repository (as shown below):

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
