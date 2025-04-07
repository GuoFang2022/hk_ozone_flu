# hk_ozone_flu

This repository contains data and codes for reproducing causal analyses on environmental determinants of influenza activity in Hong Kong via an integrative methodological framework of three disparate but complementary approaches, namely convergent cross-mapping (CCM) from a dynamical systems perspective, Peter-Clark momentary conditional independence (PCMCI) under the graphical modeling framework, and generalized linear regression model (GLM).

- 'dat_HK_fluP.rda' is fed to 'CCM.R' and 'GLM.R' R scripts to generate the main weekly results of CCM and GLM analyses, repsectively.

- 'df_pcmci_norm.csv' (normalized version of 'dat_HK_fluP.rda' in the .csv format) is fed to 'pcmci.py' Python script to run PCMCI analysis, outputing causal network graph given certain assumptions.

Note that the daily results shown in the paper during sensitivity analyses are generated using the same codes archived here, except for modifications to the lag specification and significance threshold as necessary.