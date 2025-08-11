# `sniceroclinic-coco` tracking snicercolinic with coco

This repository contains computational scripts to reproduce the a bifurcation diagram for SNICeroclinics in a 2d ODE.


## Subfolder:
* `coco_r3316`: snapshot of coco used for the runs. For newest stable version: (https://sourceforge.net/projects/cocotools/)

## Instructions

Computations require Matlab

1. Enter folder `coco_r3316`.
2. Execute `startup` script to set paths for `coco`.
3. Execute script `run_hill_js` to perform all computations, before continuation of homclinic orbit is possible. Results will be stored in subfolder `data`.
4. Execute `run_hom`. This tracks a large-period periodic orbit in two parameters. The branch spends a long time near rthe SNICeroclinic. Along the branch in the bifurcation diagram an approximate location of the saddle-node and of the saddle gets recorded.
5. Execute script `run_homsnic`, which tracks the SNICeroclinic in 3 parameters mu, gamma, beta. Results will be stored in subfolder `data/homsnic_phase` (run is called `homsnic`).
6. Execute script `plot_homsnic` to generate plots.
7. Execute script `run_ncsnic` to obtain bifurcation curve for non-central SNIC in 3 parameters mu, gamma, beta. Results will be stored in subfolder `data/ncsnic` (run is called `ncsnic`).
8. Execute script `plot_ncsnic` to generate plots.
