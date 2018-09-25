#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 15:36:33 2018

@author: carly
"""

    
#    def master_equation_handling(self, exp_dict_list:list, 
#                                 parsed_yaml_file_list:list,
#                                 master_equation_sensitivites:dict,
#                                 master_equation_reactions:list): 
#        
#        
#        def slicing_out_reactions(reaction_string,array):
#            reactions_in_cti_file = exp_dict_list[0]['simulation'].processor.solution.reaction_equations()
#            index_of_reaction_in_cti = reactions_in_cti_file.index(reaction_string)
#            column_of_array = array[:,index_of_reaction_in_cti]
#            column_of_array = column_of_array.reshape((column_of_array.shape[0],
#                                                          1))                  
#            return column_of_array
#        def assemble_slices_into_array(list_of_columns):
#            array = np.hstack((list_of_columns))
#            return array
#        def multiply_by_sensitivites(array,list_of_sensitivites):
#            sensitivity_multiplied_list = []
#            for sensitivity in list_of_sensitivites:
#                array = np.multiply(array,sensitivity)
#                array =  array.reshape((array.shape[0],1))
#                sensitivity_multiplied_list.append(array)
#            sensitivity_mapped_array = np.hstack((sensitivity_multiplied_list))
#            
#            return sensitivity_mapped_array
#        
#        def sum_of_list(list_to_sum):
#            summation = sum(list_to_sum)
#            return summation
#        def mapping_N(array,temperature):
#            for x,column in enumerate(array.T):
#                array[:,x]= np.multiply(1/np.log(temperature),column)
#            return array
#        
#        def mapping_Ea(array,temperature):
#            for x,column in enumerate(array.T):
#                array[:,x]= np.multiply(-temperature,column) * (1/8314.4621)
#            return array
#        
#        k_mapping = []   
#        for i,exp in enumerate(exp_dict_list):
#            
#            vertically_stacked_single_experiment_A = []
#            vertically_stacked_single_experiment_N = []
#            vertically_stacked_single_experiment_Ea = []
#            if parsed_yaml_file_list[i]['moleFractionObservables'][0] != None or parsed_yaml_file_list[i]['concentrationObservables'][0] !=None:
#                As = exp['ksens']['A']
#                A_observable_arrays = []
#                for observable in As:
#                    temp = []
#                    for reaction in master_equation_reactions:
#                        column = slicing_out_reactions(reaction,observable)
#                        array_mapped_to_sensitivity = multiply_by_sensitivites(column,master_equation_sensitivites[reaction]['A'])                        
#                        temp.append(array_mapped_to_sensitivity)
#                    temp = np.hstack((temp))
#
#                    A_observable_arrays.append(temp)
#                    
#                Ns = exp['ksens']['N']
#                N_observable_arrays = []    
#                for p,observable in enumerate(Ns):
#                    temp=[]
#                    for reaction in master_equation_reactions:
#                        column = slicing_out_reactions(reaction,observable)
#                        array_mapped_to_sensitivity = multiply_by_sensitivites(column,master_equation_sensitivites[reaction]['N'])
#                        #reverse mapping array
#                        array_mapped_to_N = mapping_N(array_mapped_to_sensitivity,exp['simulation'].pressureAndTemperatureToExperiment[p]['temperature'])
#                        temp.append(array_mapped_to_N)
#                    temp = np.hstack((temp))
#                    
#                    N_observable_arrays.append(temp)
#                 
#                Eas = exp['ksens']['Ea']  
#                Ea_observable_arrays = [] 
#                for observables in Eas:
#                    temp = []
#                    for reaction in master_equation_reactions:
#                        column = slicing_out_reactions(reaction,observable)
#                        array_mapped_to_sensitivity = multiply_by_sensitivites(column,master_equation_sensitivites[reaction]['Ea'])
#                        #reverse mapping of array
#                        array_mapped_to_Ea = mapping_Ea(array_mapped_to_sensitivity,exp['simulation'].pressureAndTemperatureToExperiment[p]['temperature'])                        
#                        temp.append(array_mapped_to_Ea)
#                        
#                    temp = np.hstack((temp))
#                    Ea_observable_arrays.append(temp)
#                    
#                #vertically stacked arrays
#                vertically_stacked_A_arrays = np.vstack((A_observable_arrays))
#                vertically_stacked_n_arrays = np.vstack((N_observable_arrays))
#                vertically_stacked_Ea_arrays = np.vstack((Ea_observable_arrays))
#                
#                
#                #appending things to overall experiment list which will include absorbance and 
#                #mole fractions and concentration observables if they apply
#                vertically_stacked_single_experiment_A.append(vertically_stacked_A_arrays)
#                vertically_stacked_single_experiment_N.append(vertically_stacked_n_arrays)
#                vertically_stacked_single_experiment_Ea.append(vertically_stacked_Ea_arrays)
#                
#            
#            if 'absorbance_observables' in list(exp.keys()):
#                
#                wavelengths = parsed_yaml_file_list[i]['absorbanceCsvWavelengths']
#                absorbance_k_sens = exp['absorbance_ksens']
#                #absorbance k sens is a dict with wavelenghts as keys and the values are lists of 3 things A,N,Ea sensitivity arrays
#                A_absorbance_observable_arrays = []
#                N_absorbance_observable_arrays = []
#                Ea_absorbance_observable_arrays = []
#                    
#                for k,wl in enumerate(wavelengths):
#                    temp = []
#                    for reaction in master_equation_reactions:
#                        column = slicing_out_reactions(reaction,absorbance_k_sens[wl][0])
#                        array_mapped_to_sensitivity = multiply_by_sensitivites(column,master_equation_sensitivites[reaction]['A']) 
#                        temp.append(array_mapped_to_sensitivity)
#                    temp = np.hstack((temp))
#                    A_absorbance_observable_arrays.append(temp)
#                    
#                for k,wl in enumerate(wavelengths):
#                    temp=[]
#                    for reaction in master_equation_reactions:
#                        column = slicing_out_reactions(reaction,absorbance_k_sens[wl][1])
#                        array_mapped_to_sensitivity = multiply_by_sensitivites(column,master_equation_sensitivites[reaction]['N']) 
#                        
#                        #reverse mapping of a single column 
#                        array_mapped_to_N = mapping_N(array_mapped_to_sensitivity,exp['time_history_interpolated_against_abs'][wl]['temperature'])
#                        temp.append(array_mapped_to_N)
#                    temp = np.hstack((temp))
#                    N_absorbance_observable_arrays.append(temp)
#                    
#                for k,wl in enumerate(wavelengths):
#                    temp = []
#                    for reaction in master_equation_reactions:
#                        column = slicing_out_reactions(reaction,absorbance_k_sens[wl][2])
#                        array_mapped_to_sensitivity = multiply_by_sensitivites(column,master_equation_sensitivites[reaction]['Ea'])
#                        #reverse mapping of a single column 
#                        array_mapped_to_Ea = mapping_Ea(array_mapped_to_sensitivity,exp['time_history_interpolated_against_abs'][wl]['temperature'])
#                        temp.append(array_mapped_to_Ea)
#                        
#                    temp = np.hstack(temp)
#                    Ea_absorbance_observable_arrays.append(temp)         
#                    
#                vertically_stacked_A_absorbance_array = np.vstack((A_absorbance_observable_arrays))
#                vertically_stacked_n_absorbance_array = np.vstack((N_absorbance_observable_arrays))
#                vertically_stacked_Ea_absorbance_array = np.vstack((Ea_absorbance_observable_arrays))
#                
#                vertically_stacked_single_experiment_A.append(vertically_stacked_A_absorbance_array)
#                vertically_stacked_single_experiment_N.append(vertically_stacked_n_absorbance_array)
#                vertically_stacked_single_experiment_Ea.append(vertically_stacked_Ea_absorbance_array)
#            
#            vertically_stacked_single_experiment_A = np.vstack((vertically_stacked_single_experiment_A))
#            vertically_stacked_single_experiment_N = np.vstack((vertically_stacked_single_experiment_N))
#            vertically_stacked_single_experiment_Ea = np.vstack((vertically_stacked_single_experiment_Ea))
#            
#            #mapping the single experiment back to K
#            k_mapping_single_experiment = vertically_stacked_single_experiment_A + vertically_stacked_single_experiment_N + vertically_stacked_single_experiment_Ea
#            k_mapping.append(k_mapping_single_experiment)
#            
#            
#        k_mapping = np.vstack((k_mapping))
#        
#        return k_mapping
#    
#    def surrogate_model_molecular_parameters(self,master_equation_sensitivites:dict,
#                                              master_equation_reactions:list,
#                                              delta_x_molecular_params_by_reaction_dict,
#                                              exp_dict_list):
#        
#        reactions_in_cti_file = exp_dict_list[0]['simulation'].processor.solution.reaction_equations()
#        number_of_reactions = len(reactions_in_cti_file)
#        As = []
#        Ns = []
#        Eas = []
#        for reaction in range(len(master_equation_reactions)):
#            tempA=[]
#            for i,sensitivty in enumerate(master_equation_sensitivites['R_'+str(reaction)]['A']):
#                tempA.append(sensitivty*delta_x_molecular_params_by_reaction_dict['R_'+str(reaction)][i])
#            
#            sum_A = sum(tempA)
#            As.append(sum_A)
#             
#            tempN = []
#            for i,sensitivty in enumerate(master_equation_sensitivites['R_'+str(reaction)]['N']):
#                tempA.append(sensitivty*delta_x_molecular_params_by_reaction_dict['R_'+str(reaction)][i])
#            sum_N = sum(tempN)
#            Ns.append(sum_N)
#
#            tempEa = []
#            for i,sensitivty in enumerate(master_equation_sensitivites['R_'+str(reaction)]['Ea']):
#                tempA.append(sensitivty*delta_x_molecular_params_by_reaction_dict['R_'+str(reaction)][i])
#            sum_Ea = sum(tempEa)
#            #not sure why we are multiplying by this 
#            Eas.append(sum_Ea*1000.0*1.987*4.186e3)             
#        
#       
#        
#       
#        AsNsEas =[[] for x in range(len(master_equation_reactions))]
#        for x in range(len(As)):
#            AsNsEas[x].append(As[x])
#            AsNsEas[x].append(Ns[x])
#            AsNsEas[x].append(Eas[x])
#        innerDict = ['A','n','Ea']   
#        l = [dict(zip(innerDict,AsNsEas[x])) for x in range(len(AsNsEas))]
#        Keys = []
#
#        for x in np.arange(number_of_reactions-len(master_equation_reactions),number_of_reactions):
#            Keys.append('r'+str(x))
#
#        MP = dict(zip(Keys,l))
# 
#        return MP
    

    
    