# `SNinteraction-2025` 

This repository contains interactive versions of figures and computational scripts to reproduce the results shown in figures of 

> P. Ashwin, C. Postlethwaite, J. Sieber (2026), *Local interaction of two systems with saddle-node bifurcations: mutualistic and mixed cases*. 

## Interactive figures 
(require javascript)
* [Figure2.html](https://peterashwin.github.io/SNinteraction-2025/Figure2.html)
* [Figure5.html](https://peterashwin.github.io/SNinteraction-2025/Figure5.html)
* [Figure2-Hopf-L1.html](https://peterashwin.github.io/SNinteraction-2025/Figure2-Hopf-L1.html) (not in the paper)

## Subfolders:
* [sniceroclinic-coco](sniceroclinic-coco): scripts for reproducing graphs in Figure 2, Figure 5, and check of criticality of Hopf bifurcation.
* [sniceroclinic-coco/coco_r3316](sniceroclinic-coco/coco_r3316): snapshot of coco used for the runs. For newest stable version: [sourceforge.net/projects/cocotools](https://sourceforge.net/projects/cocotools/)
* [MATLAB](MATLAB) Code to generate escape time and departure angle figures 6 and 7.  Note that plots are created for all labelled points in figures 1 and 2; only a subset of the points from figure 2 are used in figure 7 in the paper.

## Instructions

Computations require Matlab

### Folder sniceroclinic-coco

1. Enter folder [sniceroclinic-coco](sniceroclinic-coco).
2. Execute script [run_hill_js.m](sniceroclinic-coco/run_hill_js.m) to perform all computations, before continuation of homclinic orbit is possible. Results will be stored in subfolder `data`.
3. Execute [run_hom.m](sniceroclinic-coco/run_hom.m). This tracks a large-period periodic orbit in two parameters. The branch spends a long time near rthe SNICeroclinic. Along the branch in the bifurcation diagram an approximate location of the saddle-node and of the saddle gets recorded.
4. Execute script [run_homsnic.m](sniceroclinic-coco/run_homsnic.m), which tracks the SNICeroclinic in 3 parameters mu, gamma, beta. Results will be stored in subfolder `data/homsnic` (run is called `homsnic`).
5. Execute script [run_ncsnic.m](sniceroclinic-coco/run_ncsnic.m) to obtain bifurcation curve for non-central SNIC in 3 parameters mu, gamma, beta. Results will be stored in subfolder `data/ncsnic` (run is called `ncsnic`).
6. Execute script [run_TB_GH.m](sniceroclinic-coco/run_TB_GH.m) to obtain bifurcation curves for Takens-Bogdanov and Bautin points in 3 parameters mu, gamma, beta. Results will be stored in subfolders `data/hill_btp_run` and `data/hill_gh_run`.
7. Execute script [animate_hom.m](sniceroclinic-coco/animate_hom.m) to view animation of time profiles of large-period periodic orbits.
8. Execute script [plot_2dbif_mu_gamma.m](sniceroclinic-coco/plot_2dbif_mu_gamma.m) to generate plots for Figure 2 and Figure2-Hopf-L1.
9. Execute script [plot_3dbif_mu_gamma_beta.m](sniceroclinic-coco/plot_3dbif_mu_gamma_beta.m) to generate plots for Figure 5.
10. Execute script [hopf_criticality.m](sniceroclinic-coco/hopf_criticality.m) to determine formula for Lyapunov coeffficient of Hopf bifurcation.

### Folder MATLAB

Run script [hilltop_scan_v11.m](MATLAB/hilltop_scan_v11.m) to generate the figures. 
