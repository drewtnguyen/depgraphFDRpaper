This repo contains reproduction code for the figures of 
the paper "Controlling the false discovery rate under a non-parametric graphical dependence model", pre-print at <https://arxiv.org/abs/2506.24126>.

First install the `depgraphFDR` package:
```{r}
devtools::install_github("drewtnguyen/depgraphFDR", ref = "08f4329")
```
(This actually installs an old version, so that the reproduction code works even if `depgraphFDR` has been updated.)

To run the commands described below, you need to navigate to 
this repo's directory in the Terminal. 

# Simulation figures (Figs. 1, 5-8, 11-18)

## Figure 1

First, run
```{bash}
sh jobs/fig1_job.sh    
```
This will populate the `rawdata` directory with files. Next, 
run 
```{bash}
Rscript mvgauss_wrangle.R 
Rscript figure1_plot.R 
```
This will populate the `plots` directory.

## Two-sided Gaussian testing (Figs. 5, 6, 11-18)

First, run
```{bash}
sh jobs/gauss_jobs.sh    
```
or alternately 
```{bash}
parallel --jobs 8 --ungroup < jobs/gauss_jobs.sh
```
This will populate the `rawdata` directory with files. Next, 
run 
```{bash}
Rscript R/mvgauss_wrangle.R 
Rscript R/mvgauss_plots.R 
```
This will populate the `plots` directory.


## Non-null clustering (Fig. 7)

Run
```{bash}
Rscript R/thomas_pp_plots.R 
```
This will populate the `plots` directory.

## FDR inflation (Fig. 8)

First, run
```{bash}
sh jobs/adv_job.sh    
```
This will populate the `rawdata` directory with files. Next, 
run 
```{bash}
Rscript R/ctrl_wrangle.R 
Rscript R/ctrl_plots.R 
```
This will populate the `plots` directory.

# GWAS figure (Fig. 9)

## Assembling the data

The `realdata` directory should have the following structure:
```
realdata/
├── brainvar_eqtls.xlsx
├── daner_PGC_SCZ52_0513a.hq2.gz
└── g1000_eur/
    ├── g1000_eur.fam
    ├── g1000_eur.bim
    └── g1000_eur.bed
```
To obtain `brainvar_eqtls.xlsx`, download the file `media-5.xlsx` 
from [this link](https://www.biorxiv.org/content/biorxiv/early/2019/03/22/585430/DC5/embed/media-5.xlsx?download=true), which could also be found by navigating the supplement of [Werling et al 2020](https://www.biorxiv.org/content/10.1101/585430v1). Rename it to `brainvar_eqtls.xlsx`. 

To obtain `daner_PGC_SCZ52_0513a.hq2.gz`, download it from 
[this link](https://figshare.com/articles/dataset/scz2014/14672163?file=28570554), 
which could also be found by navigating the [download links of the Psychiatric Genomics Consortium.](https://pgc.unc.edu/for-researchers/download-results/)

To obtain `g1000_eur`, download the file `g1000_eur.zip` from 
[this link](https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe), 
which could also be found by navigating the [download links of the MAGMA software maintainers](https://cncr.nl/research/magma/) to find the download link corresponding 
to reference data files from Phase 3 of 1,000 Genomes (European ancestry). 


## Generate the plot

Run 
```{bash}
Rscript R/gwas_dep_plot.R 
```

