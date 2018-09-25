import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import MSI.optimization.matrix_loader as ml
import MSI.optimization.opt_runner as opt
import MSI.simulations.absorbance.curve_superimpose as csp
import MSI.simulations.yaml_parser as yp
import cantera as ct
import numpy as np
#

class MSI_shocktube_optimization(object):
        
    def __init__(self, cti_file:str,perturbment:int,
                 kineticSens:int,physicalSens:int,
                 data_directory:str,yaml_file_list:list,
                 reaction_uncertainty_csv:str,
                 k_target_values_csv:str):
        
        self.cti_file_name = cti_file
        self.perturbment = perturbment
        self.kineticSens = kineticSens
        self.physicalSens = physicalSens
        self.data_directory = data_directory
        self.yaml_file_list = yaml_file_list
        self.yaml_file_list_with_working_directory = None
        self.processor = None
        self.list_of_yaml_objects = None
        self.list_of_parsed_yamls = None
        self.experiment_dictonaries = None
        self.reaction_uncertainty_csv = reaction_uncertainty_csv
        self.k_target_values_csv = k_target_values_csv
    # call all of leis functions where we do the molecular paramter stuff and turn a flag on 
    
    def append_working_directory(self):
        yaml_file_list_with_working_directory = []
        for i, file_set in enumerate(self.yaml_file_list):
            temp = []
            for j,file in enumerate(self.yaml_file_list[i]):
                temp.append(self.data_directory+'/'+file)
            temp = tuple(temp)
            yaml_file_list_with_working_directory.append(temp)
        self.yaml_file_list_with_working_directory = yaml_file_list_with_working_directory
        return
# pre process the cti file to remove reactions and rename it,also save it as the first run of the file
        
    
    def establish_processor(self):
        processor = pr.Processor(self.data_directory +'/'+ self.cti_file_name)
        self.processor = processor
        return 
    
    def parsing_yaml_files(self,loop_counter=0):
        yaml_instance = yp.Parser()
        self.yaml_instance = yaml_instance
        list_of_yaml_objects = yaml_instance.load_yaml_list(yaml_list=self.yaml_file_list_with_working_directory)
        self.list_of_yaml_objects = list_of_yaml_objects
        list_of_parsed_yamls = yaml_instance.parsing_multiple_dictonaries(list_of_yaml_objects = list_of_yaml_objects,loop_counter=loop_counter)
        self.list_of_parsed_yamls = list_of_parsed_yamls
        #print(yaml_instance.original_experimental_conditions)
        return
    
    def running_shock_tube_simulations(self):
        optimization_instance = opt.Optimization_Utility()
        experiment_dictonaries = optimization_instance.looping_over_parsed_yaml_files(self.list_of_parsed_yamls,
                                          self.yaml_file_list_with_working_directory ,
                                          processor=self.processor, 
                                          kineticSens=self.kineticSens,
                                          physicalSens=self.physicalSens,
                                          dk=self.perturbment)
        
        
        self.experiment_dictonaries = experiment_dictonaries
        return
    
    def building_matrices(self,loop_counter=0):
        matrix_builder_instance = ml.OptMatrix()
        self.matrix_builder_instance = matrix_builder_instance
        S_matrix = matrix_builder_instance.load_S(self.experiment_dictonaries,self.list_of_parsed_yamls,dk=self.perturbment)
        self.S_matrix = S_matrix
        Y_matrix,Y_data_frame = matrix_builder_instance.load_Y(self.experiment_dictonaries,self.list_of_parsed_yamls,loop_counter=loop_counter)
        self.Y_matrix = Y_matrix
        self.Y_data_frame = Y_data_frame
        z_matrix,z_data_frame,sigma = matrix_builder_instance.build_Z(self.experiment_dictonaries,self.list_of_parsed_yamls,
                                                       loop_counter=loop_counter,
                                                       reaction_uncertainty = self.data_directory +'/'+self.reaction_uncertainty_csv )
        self.z_matrix = z_matrix
        self.z_data_frame = z_data_frame
        self.sigma = sigma
        return
    
    def adding_k_target_values(self):
        target_value_instance = ml.Adding_Target_Values(self.S_matrix,self.Y_matrix,self.z_matrix,self.sigma)
        k_target_values_for_s = target_value_instance.target_values_for_S(self.data_directory +'/'+ self.k_target_values_csv,self.experiment_dictonaries)
        k_targets_for_y = target_value_instance.target_values_Y(self.data_directory +'/'+ self.k_target_values_csv ,self.experiment_dictonaries)
        k_targets_for_z,sigma = target_value_instance.target_values_for_Z(self.data_directory +'/'+ self.k_target_values_csv)
        S_matrix,Y_matrix,z_matrix,sigma = target_value_instance.appending_target_values(k_targets_for_z,k_targets_for_y,k_target_values_for_s,sigma)
        self.S_matrix = S_matrix
        self.Y_matrix = Y_matrix
        self.z_matrix = z_matrix
        self.sigma = sigma
        return
    
    def matrix_math(self,loop_counter = 0):
        X = self.matrix_builder_instance.matrix_manipulation(loop_counter,XLastItteration = np.array(()))
        
        self.X = X
        
        
        deltaXAsNsEas,physical_observables,absorbance_coef_update_dict, X_to_subtract_from_Y = self.matrix_builder_instance.breakup_delta_x(self.X,
                                                                                                                                          self.experiment_dictonaries,
                                                                                                                                            loop_counter=loop_counter)
        self.physical_obervable_updates_list = physical_observables 
        self.absorbance_coef_update_dict = absorbance_coef_update_dict
        return
    
    def updating_files(self,loop_counter):
        updated_file_name_list = self.yaml_instance.yaml_file_updates(self.yaml_file_list_with_working_directory,
                                             self.list_of_parsed_yamls,self.experiment_dictonaries,
                                             self.physical_obervable_updates_list,
                                             loop_counter = loop_counter)
        self.updated_file_name_list = updated_file_name_list

        
        updated_absorption_file_name_list = self.yaml_instance.absorption_file_updates(self.updated_file_name_list,
                                                                                       self.list_of_parsed_yamls,
                                                                                      self.experiment_dictonaries,
                                                                                       self.absorbance_coef_update_dict,
                                                                                       loop_counter = loop_counter)
        self.updated_absorption_file_name_list = updated_absorption_file_name_list
        
        #update the cti files pass in the renamed file 
        #self.processor.cti_write2(x={},original_cti='',master_rxns='',master_index=[],MP={})
        
        return
    
    def one_run_shock_tube_optimization(self,loop_counter=0):
        self.append_working_directory()
        self.establish_processor()
        self.parsing_yaml_files(loop_counter = loop_counter)
        if loop_counter == 0:
            original_experimental_conditions_local = self.yaml_instance.original_experimental_conditions
            self.original_experimental_conditions_local = original_experimental_conditions_local
        self.running_shock_tube_simulations()
        self.building_matrices(loop_counter=loop_counter)
        self.adding_k_target_values()
        self.matrix_math(loop_counter=loop_counter)
        self.updating_files(loop_counter=loop_counter)
        
        
        
    def multiple_shock_tube_runs(self,loops):
        for loop in range(loops):
            self.one_run_shock_tube_optimization(loop_counter=loop)
        
        return
            
    
                                                                  
        
                                                       
            
                



