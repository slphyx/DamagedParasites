## Overview

This repository contains the code used in the publication:

**An artesunate pharmacometric model to explain therapeutic responses in *falciparum* malaria**  
Saralamba S. *et al.*  
*Journal of Antimicrobial Chemotherapy*, Volume 78, Issue 9, September 2023, Pages 2192–2202  
[https://doi.org/10.1093/jac/dkad219](https://doi.org/10.1093/jac/dkad219)

## Citation

If you use this code in your work, please cite the original publication:

Saralamba S, et al. *An artesunate pharmacometric model to explain therapeutic responses in falciparum malaria.*  
Journal of Antimicrobial Chemotherapy. 2023;78(9):2192–2202.  
doi:[10.1093/jac/dkad219](https://doi.org/10.1093/jac/dkad219)

You may also cite this repository directly:

@misc{DamagedParasites, 
author = {Saralamba, S. and collaborators}, 
title = {DamagedParasites: Code for artemisinin resistance mathematical model}, 
year = {2023}, 
publisher = {GitHub}, 
journal = {GitHub repository}, 
howpublished = {\url{https://github.com/slphyx/DamagedParasites}} }

## Files	
	|- C# -
	|     |- DamagedParasites_C#.zip  => the proposed model written in C#
	|
	|- Stan - 
	|	|- recoveryPailin.stan	=> the proposed model written in Stan and the code for fitting the model to the data
	|
	|- Wolfram -
	           |- Figure2and4-PL.nb   =>  code for generating Figure 2 and 4
		   |- Figure2and4-WP.nb	  =>  code for generating Figure 2 and 4
		   |- Figure3.nb          =>  code for generating Figure 3
		   |- Figure5.nb          =>  code for generating Figure 5
		   |- Figure6.nb          =>  code for generating Figure 6 
		   |- PL-stan-output.zip  =>  Stan output for Pailin
		   |- WP-stan-output.zip  =>  Stan output for Wang Pha  
		   |- PLsummary.csv       =>  Stan output summary
		   |- WPsummary.csv       =>  Stan output summary
			   
