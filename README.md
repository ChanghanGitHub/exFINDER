# exFINDER: Identifying external communication signals using single-cell transcriptomics data

## Introduction
exFINDER is a method that identifies external signals (received signals from the external system) in the single-cell transcriptomics datasets by utilizing the prior knowledge of signaling pathways. Specifically,exFINDER contains the following features:

- It develops a computational method that links the prior knowledge and scRNA-seq data to identify the external signals received by the cells. 
- It uses a computational model that utilizes the graph theory to infer the external signal-target signaling networks (exSigNet) and link the expression data to signaling strength. 
- It provides multiple ways for visualization and analysis of the external signal-target signaling networks (exSigNet).


## Installation

exFINDER R package can be easily installed from Github using devtools (it may take a little while):  

```
devtools::install_github("ChanghanGitHub/exFINDER")
```

## Applications of exFINDER on different datasets
- [exFINDER identifies differentiation-associated external signals during zebrafish neural crest (NC) development. ](https://htmlpreview.github.io/?https://github.com/ChanghanGitHub/exFINDER/blob/master/vignettes/vignette1_Tatarakis2021_Zebrafish.html), source data comes from [Tatarakis et al., Cell reports, 2021](https://www.sciencedirect.com/science/article/pii/S2211124721016363)

- [exFINDER suggests critical external signals and targets during sensory neurogenesis in mouse. ](https://htmlpreview.github.io/?https://github.com/ChanghanGitHub/exFINDER/blob/master/vignettes/vignette2_Faure2020_Mouse_Part1.html), source data comes from [Faure et al., Nature communications, 2020](https://www.nature.com/articles/s41467-020-17929-4)

- [exFINDER predicts the roles of external signals and uncovers transition paths in differentiation](https://htmlpreview.github.io/?https://github.com/ChanghanGitHub/exFINDER/blob/master/vignettes/vignette3_Faure2020_Mouse_Part2.html), source data comes from [Faure et al., Nature communications, 2020](https://www.nature.com/articles/s41467-020-17929-4)


## Benchmarking
- [Benchmarking exFINDER using human skin data.](https://htmlpreview.github.io/?https://github.com/ChanghanGitHub/exFINDER/blob/master/vignettes/vignette0_Benchmarking_He2020_Human.html), source data comes from [He, Helen, et al. Journal of Allergy and Clinical Immunology, 2020](https://www.sciencedirect.com/science/article/pii/S0091674920301822?casa_token=CsRgPXXj_M8AAAAA:MldVWJ0Oug2g-aXSg86-De9IJcq3TZ0v9TyhQT2kDWwZa4msBLR4oFsU7LhvJ_ZQsGrwOY_pbCM)




