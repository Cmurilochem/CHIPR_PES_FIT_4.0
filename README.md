# CHIPR-4.0

CHIPR-4.0 is a general software to fit global potential energy surfaces of diatomic, triatomic and tetratomic molecules using ab initio data points as calibrating set and the Combined-Hyperbolic-Inverse-Power-Representation (CHIPR) method.

For diatomic molecules, the code also allows users obtain experimentally-derived potential energy curves by performing a direct-fit to available spectroscopic data. 

The program performs an automatic global minimum search on the final fitted PES, harmonic vibrational analysis on the global minimum, in addition to provide ready-to-use subroutines containing the fitted n-body terms. These are output as CHIPR_DIAT_FUNC.f90, CHIPR_TRIAT_FUNC.f90 and CHIPR_TETRA_FUNC.f90 for two-body, three-body and four-body terms, respectively.

The source code can be found in "CHIPR-4.0_SOURCE_CODE". Test runs for C<sub>2</sub>, N<sub>2</sub>, SiC, C<sub>3</sub>, SiC<sub>2</sub>, C<sub>3</sub>H, and C<sub>2</sub>H<sub>2</sub> are provided as separate directories. 

see https://data.mendeley.com/datasets/8wdv87gt5x/2 to access the original open source version of the code and publication.

Authors:
Carlos M. R. Rocha* (carlosmurilorocha@gmail.com) and Antonio J. C. Varandas (varandas@uc.pt)

-----------------------------------------------------------------------------------------------------------------------------------------------------------

When using this program, check the following references:

(1). For the CHIPR method:

   A. J. C. Varandas, J. Chem. Phys. 138, 054120 (2013); https://doi.org/10.1063/1.4788912
 
   A. J. C. Varandas, J. Chem. Phys. 138, 134117 (2013); https://doi.org/10.1063/1.4795826
 
   C. M. R. Rocha and A. J. C. Varandas, J. Phys. Chem. A 123, 38 (2019); https://doi.org/10.1021/acs.jpca.9b03194
 
   C. M. R. Rocha and A. J. C. Varandas, Phys. Chem. Chem. Phys. 21, 24406 (2019); https://doi.org/10.1039/c9cp04890a
 
   C. M. R. Rocha, H. Linnartz and A. J. C. Varandas, J. Chem. Phys. 157, 104301 (2022); https://doi.org/10.1063/5.0096364
 
(2). For the code: 

   C. M. R. Rocha and A. J. C. Varandas, Comput. Phys. Commun. 258, 107556 (2021); https://doi.org/10.1016/j.cpc.2020.107556

