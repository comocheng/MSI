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
import MSI.utilities.plotting_script as plotter

#
#MSI_st_instance_one = stMSI.MSI_shocktube_optimization('FFCM1_custom_test.cti',
#                                                   .01,
#                                                   1,
#                                                   1,
#                                                   'MSI/data/test_data',
#                                                   [['Hong_1.yaml'],
#                                                    ['Hong_2.yaml'],
#                                                    ['Hong_3.yaml'],
#                                                    ['Hong_4.yaml','Hong_4_abs.yaml'],
#                                                    ['Troe_4.yaml','Troe_4_abs.yaml'],
#                                                    ['Troe_5.yaml','Troe_5_abs.yaml'],
#                                                    ['Troe_6.yaml','Troe_6_abs.yaml'],
#                                                    ['Troe_7.yaml','Troe_7_abs.yaml'],
#                                                    ['Troe_8.yaml','Troe_8_abs.yaml'],
#                                                    ['Hong_5.yaml','Hong_5_abs.yaml']],                 
#                                                   'uncertainty_test.csv','FFCM1_custom_target_value_test.csv' )
                 
                 
                 

MSI_st_instance_one = stMSI.MSI_shocktube_optimization('FFCM1_custom_original.cti',
                                                   .01,
                                                   1,
                                                   1,
                                                   'MSI/data/test_data',
                                                   [['Hong_1.yaml']],                 
                                                   'uncertainty_test.csv','' )
MSI_st_instance_one.one_run_shock_tube_optimization()

#Ydf = MSI_st_instance_one.Y_data_frame
#S_matrix_original = MSI_st_instance_one.S_matrix
#exp_dict_list_original = MSI_st_instance_one.experiment_dictonaries
#original_covariance = MSI_st_instance_one.original_covarience
#X_to_subtract = MSI_st_instance_one.X_to_subtract_from_Y
#
S_matrix = MSI_st_instance_one.S_matrix
#s_matrix = MSI_st_instance_one.S_matrix
##need to fix this and return _s_matrix and y_matrix
#Y_matrix = MSI_st_instance_one.Y_matrix
##y_matrix = MSI_st_instance_one.y_matrix
#z_matrix = MSI_st_instance_one.z_matrix
#Zdf = MSI_st_instance_one.z_data_frame
#X = MSI_st_instance_one.X
#covarience = MSI_st_instance_one.covarience
#
#parsed_yaml_list = MSI_st_instance_one.list_of_parsed_yamls
#sigma =MSI_st_instance_one.sigma


exp_dict_list_optimized = MSI_st_instance_one.experiment_dictonaries
timeHistory = exp_dict_list_optimized[0]['simulation'].timeHistories[0]
experimental_data = exp_dict_list_optimized[0]['experimental_data']
physical_obervable_updates_list = MSI_st_instance_one.physical_obervable_updates_list
As_nS_eas_updates = MSI_st_instance_one.deltaXAsNsEas
#deltaXAsNsEas = MSI_st_instance_one.deltaXAsNsEas


#MSI_st_instance_two = stMSI.MSI_shocktube_optimization('FFCM1_custom_original.cti',
#                                                   .01,
#                                                   1,
#                                                   1,
#                                                   'MSI/data/test_data',
#                                                   [['Hong_2.yaml']],                 
#                                                   'uncertainty_test.csv','' )
#MSI_st_instance_two.multiple_shock_tube_runs(2)
#S_matrix = MSI_st_instance_two.S_matrix
#s_matrix = MSI_st_instance_two.s_matrix
#Y_matrix = MSI_st_instance_two.Y_matrix
#y_matrix = MSI_st_instance_two.y_matrix
#z_matrix = MSI_st_instance_two.z_matrix
#X = MSI_st_instance_two.X
#covarience = MSI_st_instance_two.covarience
#exp_dict_list_optimized = MSI_st_instance_two.experiment_dictonaries
#parsed_yaml_list = MSI_st_instance_two.list_of_parsed_yamls
#target_value_rate_constant_csv = 'MSI/data/test_data/FFCM1_custom_target_value_test.csv'









#plotting_instance = plotter.Plotting(S_matrix,
#                                     s_matrix,
#                                     Y_matrix,
#                                     Y_matrix,
#                                     z_matrix,
#                                     X,
#                                     sigma,
#                                     covarience,
#                                     original_covariance,
#                                     S_matrix_original,
#                                     exp_dict_list_optimized,
#                                     exp_dict_list_original,
#                                     parsed_yaml_list,
#                                     '')
#plotting_instance.plotting_observables()
