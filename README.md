# Scripts that generate the results in SAASI

## Required packages/software:
R (Version 4.4.1)
Python (Version 3.13.7)

TREETIME(Python)
PASTML(Python)

## R library used
tidyr
stringr
treeio
tidyverse
ggalluvial
ggpubr
irr
stats4
ape
deSolve
strap
ips
ggplot2
ggtree
reshape2
ggstance
diversitree
ggimage
patchwork
janitor
phytools
dplyr
rstatix
jsonlite
rootSolve

## Content
### Code:
Code for generating the figures in the manuscript.

### Data:
H5N1 tree used in the analysis. See original paper (Nguyen et al. 2025, DOI: 10.1126/science.adq0900
) and GitHub: https://github.com/flu-crew/dairy-cattle-hpai-2024/tree/main/treetime.

### Figures:
Figures in the manuscript.

### SAASI:
The algorithm that is described in the manuscript.

## Running Shell script

To run the shell script, run:
```
./script_1.sh 1000 5000 5000 folder
```
The first argument is the number of simulations.
The second argument is the maximum running time.
The third argument is the maximum number of tips in the present day.
The forth argument is the folder's name


