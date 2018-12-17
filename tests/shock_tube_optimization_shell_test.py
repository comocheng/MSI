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
 #burke_target_value_test.csv                
files_to_include = [['Hong_single_data_point.yaml']]
                    
              
numer_of_iterations = 3
     

MSI_st_instance_one = stMSI.MSI_shocktube_optimization('chem_original_burke.cti',
                                                   .01,
                                                   1,
                                                   1,
                                                   'MSI/data/test_data',
                                                   files_to_include,                 
                                                   'burke_uncertainty_test.csv','' )
MSI_st_instance_one.one_run_shock_tube_optimization()

S_matrix_original = MSI_st_instance_one.S_matrix
exp_dict_list_original = MSI_st_instance_one.experiment_dictonaries
original_covariance = MSI_st_instance_one.covarience
X_one_itteration = MSI_st_instance_one.X
#need to fix this and return _s_matrix and y_matrix



MSI_st_instance_two = stMSI.MSI_shocktube_optimization('chem_original_burke.cti',
                                                   .01,
                                                   1,
                                                   1,
                                                   'MSI/data/test_data',
                                                   files_to_include,                 
                                                   'burke_uncertainty_test.csv','' )
X_list = MSI_st_instance_two.multiple_shock_tube_runs(numer_of_iterations)


deltaXAsNsEas = MSI_st_instance_two.deltaXAsNsEas
physical_obervable_updates_list = MSI_st_instance_two.physical_obervable_updates_list
absorbance_observables_updates_list = MSI_st_instance_two.absorbance_coef_update_dict
Ydf = MSI_st_instance_two.Y_data_frame
Zdf = MSI_st_instance_two.z_data_frame
experimental_dicts = MSI_st_instance_two.experiment_dictonaries
z_matrix = MSI_st_instance_two.z_matrix
s_matrix = MSI_st_instance_two.s_matrix
y = MSI_st_instance_two.y_matrix
Y_matrix = MSI_st_instance_two.Y_matrix
S_matrix = MSI_st_instance_two.S_matrix

X = MSI_st_instance_two.X
covarience = MSI_st_instance_two.covarience
exp_dict_list_optimized = MSI_st_instance_two.experiment_dictonaries
parsed_yaml_list = MSI_st_instance_two.list_of_parsed_yamls
sigma = MSI_st_instance_two.sigma
X = MSI_st_instance_two.X
delta_X = MSI_st_instance_two.delta_X
#target_value_rate_constant_csv = 'MSI/data/test_data/FFCM1_custom_target_value_test.csv'
original_cti_file = MSI_st_instance_two.data_directory +'/'+ MSI_st_instance_two.cti_file_name

experiment_dict_uncertainty = MSI_st_instance_two.experiment_dict_uncertainty_original
target_value_csv = MSI_st_instance_two.data_directory +'/'+ MSI_st_instance_two.k_target_values_csv


#k_target_value_S_matrix = MSI_st_instance_two.k_target_values_for_s


##########################################################################################################################
#PLOTTING##
##########################################################################################################################


plotting_instance = plotter.Plotting(S_matrix,
                                     s_matrix,
                                     Y_matrix,
                                     Y_matrix,
                                     z_matrix,
                                     X,
                                     sigma,
                                     covarience,
                                     original_covariance,
                                     S_matrix_original,
                                     exp_dict_list_optimized,
                                     exp_dict_list_original,
                                     parsed_yaml_list,
                                     Ydf,
                                     target_value_rate_constant_csv= MSI_st_instance_two.data_directory +'/'+'burke_target_value_single_reactions.csv' ,
                                     k_target_value_S_matrix =None,
                                     k_target_values='Off')

observable_counter_and_absorbance_wl,length_of_experimental_data = plotting_instance.lengths_of_experimental_data()
sigmas_optimized,test = plotting_instance.calculating_sigmas(S_matrix,covarience)
sigmas_original,test2 = plotting_instance.calculating_sigmas(S_matrix_original,original_covariance)
plotting_instance.plotting_observables(sigmas_original = sigmas_original,sigmas_optimized= sigmas_optimized)
diag = plotting_instance.getting_matrix_diag(covarience)
plotting_instance.Y_matrix_plotter(Y_matrix,exp_dict_list_optimized,y,sigma)



plotting_instance.plotting_rate_constants(optimized_cti_file=MSI_st_instance_two.new_cti_file,
                                original_cti_file=original_cti_file,
                                initial_temperature=250,
                                final_temperature=2500)
                                


sensitivity, top_sensitivity = plotting_instance.sort_top_uncertainty_weighted_sens()
obs = plotting_instance.plotting_uncertainty_weighted_sens()


#plotting_instance.plotting_X_itterations(list_of_X_values_to_plot = [0,1,2,3,4,5,50],list_of_X_array=X_list,number_of_iterations=numer_of_iterations)

