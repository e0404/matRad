---------------------------------------------------
Neutron RBE Data
---------------------------------------------------
A)  Constant factor for neutron RBE in file 'neutronRBE_constFactor.txt' is 
    loaded using the plan variable: 
    pln.propOpt.bioOptimization = 'const_RBExD_n'.
    Ref.: The value is set to 3 taken from Specht et al. (2015) p. 2.
---------------------------------------------------
B)  Tabulated values are taken from ICRP 103 and should actually not be 
    used for therapeutic dose calculations but for radiation protection
    purposes. With this weighting factors the equivalent dose is calcuated.
    In tabulated form, RBE values have to be given for the same energy
    intervals like the used KERMA values. Otherwise result will be wrong.
    Variable neutron RBE in file 'neutronRBE_ICRP103.txt' is 
    loaded using the plan variable: 
    pln.propOpt.bioOptimization = 'var_RBExD_n_ICRP103'.
    Ref.: ICRP Publication 103
---------------------------------------------------
C)  Tabulated values for RBE are taken from plot in Steward et al. 2015. 
    Values are interpolated in order to match KERMA dose intervals given 
    by predefined KERMA values. Note that reasonable cut-off energies
    were defined s.th. tabulated values below cut-off were averaged.