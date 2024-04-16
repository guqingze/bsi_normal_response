# Tracking normal BSI response

R code for paper "Distinct patterns of vital sign and inflammatory marker responses in adults with suspected bloodstream infection". 
DOI:https://doi.org/10.1016/j.jinf.2024.106156

## Input data
The data analysed are available from the Infections in Oxfordshire Research Database (https://oxfordbrc.nihr.ac.uk/research-themes/modernising-medical-microbiology-and-big-infection-diagnostics/iord-about/), subject to an application and research proposal meeting on the ethical and governance requirements of the Database.

## Software and packages
Analyses were performed using statistical software R, version 4.1.0 (R Project for Statistical Computing). Model fitting was performed using the ‘lme4’ package (version 1.1-27.1) (using restricted maximum likelihood estimation (REML) with the ‘t-tests use Satterthwaite’s method’ approach for approximating degrees of freedom, as implemented in the ‘lmerModLmerTest’ function), the ‘lcmm’ package (version 2.1.0), and the ‘gamlss’ package (version 5.4-20).

## Code
- linear_mixed_models: R code for linear mixed model analysis, examining the association clinical response trajectories and both infection sources and different pathogen groups.
- latent_class_mixed_models: R code for latent class mixed model analysis, identifying underlying heterogeneity in CRP response trajectories.
- lms: R code for generalized additive mixed model analysis, developing CRP centile charts characterising typical responses to guide individualised management. 