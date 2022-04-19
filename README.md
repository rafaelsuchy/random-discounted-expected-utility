# Supplementary Material: </br> Random Models for the Joint Treatment of Risk and Time Preferences

This repository contains the supplementary material for the paper "Random Models for the Joint Treatment of Risk and Time Preferences" (2022) by Jose Apesteguia, Miguel A. Ballester and Ángelo Gutierrez-Daza.

### Overview

We use the original data from [Andersen et al. (2008, AHLR)](https://www.jstor.org/stable/40056458) ([download link](https://www.econometricsociety.org/publications/econometrica/2008/05/01/eliciting-risk-and-time-preferences)) and [Andreoni and Sprenger (2012, AS)](https://www.aeaweb.org/articles?id=10.1257/aer.102.7.3357) ([download link](https://www.openicpsr.org/openicpsr/project/112570/version/V1/view)). Once the data is downloaded, it should be stored as .xls and .dta-files in the folder [`Data/raw_data`](/Data/raw_data/) under the following names: 
- `AHLR_raw_data.dta` and `AHLR_raw_data.xls`
- `AS_raw_data.dta` and `AS_raw_data.xls`. 

The folder [`Codes`](/Codes) contains the programs to replicate the empirical results that are reported in the paper. 

### Replication 

***Step 1:*** Prepare the [`Data/raw_data`](/Data/raw_data/) for the estimation routines. This can be done by running the programs [**AHLR_preparation_data.R**](/Codes/AHLR_data_preparation.R) and [**AS_data_preparation.R**](/Codes/AS_data_preparation.R), respectively. The output is stored in the [`Data/processed_data`](Data/processed_data). The processed data is used by the estimation programs. For convenience, the prepared data is already stored in [`Data/processed_data/AHLR`](Data/processed_data/AHLR) and [`Data/processed_data/AS`](Data/processed_data/AS).

***Step 2:*** Estimate the Multiple Price Lists and Convex Menus models:
- Table 1, Table 3, Table 7, Table 9, Figure 1, and Figure 2 are produced by [**AHLR_estimation.R**](/Codes/AHLR_estimation.R)
- Table 2, Table 4 and Figure 3 are produced by [**AS_estimation.R**](/Codes/AS_estimation.R)
- Table 5 is produced by [**AHLR_estimation_tremble.R**](/Codes/AHLR_estimation_tremble.R)
- Table 6 is produced by [**AS_estimation_tremble.R**](/Codes/AS_estimation_tremble.R)
- Table 8 is produced by [**AHLR_estimation_hyperbolic.R**](/Codes/AHLR_estimation_hyperbolic.R)

The estimation output for the tables is stored in the folder [`Output/estimates`](/Output/estimates). Figures are stored in the [`Output/figures`](/Output/figures) folder. 

### Main References

- [Andersen, S., Harrison, G. W., Lau, M. I., & Rutström, E. E. (2008). Eliciting Risk and Time Preferences. *Econometrica, 76*(3), 583–618.](https://www.jstor.org/stable/40056458?seq=1)
- [Andreoni, J., Sprenger, C. (2012). Risk Preferences Are Not Time Preferences. *American Economic Review, 102*(7), 3357-76.](https://www.aeaweb.org/articles?id=10.1257/aer.102.7.3357)

