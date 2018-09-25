import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import MSI.optimization.matrix_loader as ml
import MSI.optimization.opt_runner as opt
import MSI.simulations.absorbance.curve_superimpose as csp
import MSI.simulations.yaml_parser as yp
import MSI.optimization.shock_tube_optimization_shell as stMSI
import cantera as ct



                 
                 
                 

MSI_st_instance = stMSI.MSI_shocktube_optimization('FFCM1.cti',.01,1,1,'MSI/data/test_data',
                                 [['Hong_4.yaml','Hong_4_abs.yaml'],['Hong_1.yaml'],
                  ['Troe_6.yaml','Troe_6_abs.yaml']],                 
                   'uncertainty_test.csv','FFCM1_target_values.csv' )
MSI_st_instance.one_run_shock_tube_optimization()
X = MSI_st_instance.X
