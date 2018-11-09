import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cantera as ct

class Plotting(object):
    def __init__(self,S_matrix,
                 s_matrix,
                 Y_matrix,
                 y_matrix,
                 z_matrix,
                 X,
                 sigma,
                 covarience,
                 original_covariance,
                 S_matrix_original,
                 exp_dict_list_optimized,
                 exp_dict_list_original,
                 parsed_yaml_list,
                 target_value_rate_constant_csv=''):
        self.S_matrix = S_matrix
        self.s_matrix = s_matrix
        self.Y_matrix = Y_matrix
        self.y_matrix = y_matrix
        self.z_matrix = z_matrix
        self.X = X
        #self.sigma = sigma
        self.covarience=covarience
        self.original_covariance=original_covariance
        #original
        self.S_matrix_original=S_matrix_original
        self.exp_dict_list_optimized = exp_dict_list_optimized
        self.exp_dict_list_original = exp_dict_list_original
        self.parsed_yaml_list = parsed_yaml_list
        self.target_value_rate_constant_csv = target_value_rate_constant_csv
        
        
        
 #fix all the indexing to have a captial or lowercase time situation or add the module that lets you do either to all the scripts  

    def lengths_of_experimental_data(self):
        simulation_lengths_of_experimental_data = []
        for i,exp in enumerate(self.exp_dict_list_optimized):
            length_of_experimental_data=[]
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                if observable in exp['mole_fraction_observables']:
                    length_of_experimental_data.append(exp['experimental_data'][observable_counter]['Time'].shape[0])
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:
                    length_of_experimental_data.append(exp['experimental_data'][observable_counter]['Time'].shape[0])
                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                absorbance_wl=0
                for k,wl in enumerate(wavelengths):
                    length_of_experimental_data.append(exp['absorbance_experimental_data'][k]['time'].shape[0])
                    absorbance_wl+=1
                    
            simulation_lengths_of_experimental_data.append(length_of_experimental_data)
            
                    
        self.simulation_lengths_of_experimental_data=simulation_lengths_of_experimental_data
        return observable_counter+absorbance_wl,length_of_experimental_data
                    

        
    def calculating_sigmas(self,simulation_lengths_of_experimental_data,S_matrix,covarience):        
        sigmas =[[] for x in range(len(simulation_lengths_of_experimental_data))]
                 
        counter=0
        for x in range(len(simulation_lengths_of_experimental_data)):
            for y in range(len(simulation_lengths_of_experimental_data[x])):
                temp=[]
                for z in np.arange(counter,(simulation_lengths_of_experimental_data[x][y]+counter)):                             
                    SC = np.dot(S_matrix[z,:],covarience)
                    sigma = np.dot(SC,np.transpose(self.S_matrix[z,:]))
                    sigma = np.sqrt(sigma)
                    temp.append(sigma)
                temp = np.array(temp)            
                sigmas[x].append(temp)
        
                
                counter = counter + simulation_lengths_of_experimental_data[x][y]
        
        return sigmas
    
    
    
    def plotting_observables(self,sigmas_original=[],sigmas_optimized=[]):
        #
        
        
        
        for i,exp in enumerate(self.exp_dict_list_optimized):
            observable_counter=0
            for j,observable in enumerate(exp['mole_fraction_observables'] + exp['concentration_observables']):
                if observable == None:
                    continue
                plt.figure()
                if observable in exp['mole_fraction_observables']:
                    plt.plot(exp['simulation'].timeHistories[0]['time']*1e3,exp['simulation'].timeHistories[0][observable],'b',label='MSI')
                    plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original['simulation'].timeHistories[0][observable],'r',label= "$\it{A priori}$ model")
                    plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,exp['experimental_data'][observable_counter][observable],'o',color='black',label='Experimental Data')
                    plt.xlabel('Time (ms)')
                    plt.ylabel('Mole Fraction'+''+str(observable))
                    plt.title('Experiment_'+str(i+1))
                    
                    
                    

                    
                    if bool(sigmas_optimized) == True:
                        
                        high_error_optimized = np.exp(sigmas_optimized[i][observable_counter])                   
                        high_error_optimized = np.multiply(high_error_optimized,exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                        low_error_optimized = np.exp(sigmas_optimized[i][observable_counter]*-1)
                        low_error_optimized = np.multiply(low_error_optimized,exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,  high_error_optimized,'b','--')
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,low_error_optimized,'b','--')
                        
                        
                        
                        high_error_original = np.exp(sigmas_original[i][observable_counter])
                        high_error_original = np.multiply(high_error_original,self.exp_dict_list_original[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                        low_error_original = np.exp(sigmas_original[i][observable_counter]*-1)
                        low_error_original = np.multiply(low_error_original,self.exp_dict_list_original[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,  high_error_original,'--','r')
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,low_error_original,'--','r')
                    
                    plt.savefig('Experiment_'+str(i+1)+'_'+str(observable)+'.pdf', bbox_inches='tight')
                    
                    observable_counter+=1
                    
                if observable in exp['concentration_observables']:
                    plt.plot(exp['simulation'].timeHistories[0]['time']*1e3,exp['simulation'].timeHistories[0][observable]*1e6,'b',label='MSI')
                    plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original[i]['simulation'].timeHistories[0][observable]*1e6,'r',label= "$\it{A priori}$ model")
                    plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,exp['experimental_data'][observable_counter][observable+'_ppm'],'o',color='black',label='Experimental Data') 
                    plt.xlabel('Time (ms)')
                    plt.ylabel('ppm'+''+str(observable))
                    plt.title('Experiment_'+str(i+1))
                    
                    if bool(sigmas_optimized)==True:
                        high_error_optimized = np.exp(sigmas_optimized[i][observable_counter])                   
                        high_error_optimized = np.multiply(high_error_optimized,exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                        low_error_optimized = np.exp(sigmas_optimized[i][observable_counter]*-1)
                        low_error_optimized = np.multiply(low_error_optimized,exp['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,  high_error_optimized,'b','--')
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,low_error_optimized,'b','--')                    
                        
    
    
                        high_error_original = np.exp(sigmas_original[i][observable_counter])
                        high_error_original = np.multiply(high_error_original,self.exp_dict_list_original[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                        low_error_original = np.exp(sigmas_original[i][observable_counter]*-1)
                        low_error_original = np.multiply(low_error_original,self.exp_dict_list_original[i]['simulation'].timeHistoryInterpToExperiment[observable].dropna().values*1e6)
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,  high_error_original,'--','r')
                        plt.plot(exp['experimental_data'][observable_counter]['Time']*1e3,low_error_original,'--','r')
                    
                    plt.savefig('Experiment_'+str(i+1)+'_'+str(observable)+'.pdf', bbox_inches='tight')
                    observable_counter+=1
                    

            if 'perturbed_coef' in exp.keys():
                wavelengths = self.parsed_yaml_list[i]['absorbanceCsvWavelengths']
                plt.figure()
                for k,wl in enumerate(wavelengths):
                    plt.plot(exp['simulation'].timeHistories[0]['time']*1e3,exp['absorbance_calculated_from_model'][wl],'b',label='MSI')
                    plt.plot(self.exp_dict_list_original[i]['simulation'].timeHistories[0]['time']*1e3,self.exp_dict_list_original[i]['absorbance_calculated_from_model'][wl],'r',label= "$\it{A priori}$ model")
                    plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,exp['absorbance_experimental_data'][k]['Absorbance_'+str(wl)],'o',color='black',label='Experimental Data')
                    plt.xlabel('Time (ms)')
                    plt.ylabel('Absorbance'+''+str(wl))
                    plt.title('Experiment_'+str(i+1))
                    
                    if bool(sigmas_optimized)==True:
                        high_error_optimized = np.exp(sigmas_optimized[i][observable_counter])
                        high_error_optimized = np.multiply(high_error_optimized,exp['absorbance_model_data'][wl])
                        low_error_optimized = np.exp(sigmas_optimized[i][observable_counter]*-1)
                        low_error_optimized = np.multiply(low_error_optimized,exp['absorbance_model_data'][wl])
                        
                        plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,high_error_optimized,'b','--')
                        plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,low_error_optimized,'b','--')
                        
                        high_error_original = np.exp(sigmas_original[i][observable_counter])
                        high_error_original = np.multiply(high_error_original,self.exp_dict_list_original['absorbance_model_data'][wl])
                        low_error_original =  np.exp(sigmas_original[i][observable_counter]*-1)
                        low_error_original = np.multiply(low_error_original,self.exp_dict_list_original['absorbance_model_data'][wl])
                        
                        plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,high_error_original,'r','--')
                        plt.plot(exp['absorbance_experimental_data'][k]['time']*1e3,low_error_original,'r','--')
                    
                    #start here
                    plt.savefig('Experiment_'+str(i+1)+' '+'Absorbance at'+'_'+str(wl)+'.pdf', bbox_inches='tight')
                    
                    
                    

# make function to plot rate constants 
                    
    def plotting_rate_constants(self,optimized_cti_file='',
                                original_cti_file='',
                                initial_temperature=250,
                                final_temperature=2500,
                                target_value_rate_constant_csv='',
                                reactions_to_plot=[]):
        
        gas_optimized = ct.Solution(optimized_cti_file)
        gas_original = ct.Solution(original_cti_file)
        def unique_list(seq):
            checked = []
            for e in seq:
                if e not in checked:
                    checked.append(e)
            return checked
        
        def sort_rate_constant_target_values(parsed_csv,unique_reactions,gas):
            reaction_list_from_mechanism = gas.reaction_equations()
            target_value_ks = [[] for reaction in len(unique_reactions)]
            target_value_temps = [[] for reaction in len(unique_reactions)]
            reaction_list_from_mechanism = gas.reaction_equations()
            
            for i,reaction in enumerate(parsed_csv['Reaction']):
                idx = reaction_list_from_mechanism.index(reaction)
                target_value_ks[unique_reactions.index(idx)].append(parsed_csv['k'][i])
                target_value_temps[unique_reactions.index(idx)].append(parsed_csv['temperature'][i])
                
            return target_value_temps,target_value_ks
        def rate_constant_over_temperature_range_from_cantera(reaction_number,
                                                              gas,
                                                              initial_temperature=250,
                                                              final_temperature=2500,
                                                              pressure=1,
                                                              conditions = {'H2':2,'O2':1,'N2':4}):
            Temp = []
            k = []
            for temperature in np.arange(initial_temperature,final_temperature,1):
                gas.TPX = temperature,pressure*101325,conditions
                Temp.append(temperature)
                k.append(gas.forward_rate_constants[reaction_number])
            
            return Temp,k

        def calculate_sigmas_for_rate_constants(s_matrix,k_target_values_parsed_csv,unique_reactions,gas,covarience):
            #grab the last x rows of the matricies that we need 
            #start here tomorrow 
            
            reaction_list_from_mechanism = gas.reaction_equations()
            sigma_list_for_target_ks = [[] for reaction in len(unique_reactions)]
            shape = k_target_values_parsed_csv.shape
            k_values_in_s_matrix = s_matrix[-shape[0]:,:]
            for row in range(shape(k_values_in_s_matrix[0])):
                SC = np.dot(k_values_in_s_matrix[row,:],covarience)
                sigma_k = np.dot(SC,np.transpose(k_values_in_s_matrix[row,:]))
                sigma_k = np.sqrt(sigma_k)
                indx = reaction_list_from_mechanism.index(k_target_values_parsed_csv['Reaction'][row])
                sigma_list_for_target_ks[unique_reactions.index(indx)].append(sigma_k)
                
            return sigma_list_for_target_ks
        
        def calculating_target_value_ks_from_cantera_for_sigmas(k_target_values_parsed_csv,gas,unique_reactions):
            target_value_ks = [[] for reaction in len(unique_reactions)]
            
            target_reactions = k_target_values_parsed_csv['Reaction']
            target_temp = k_target_values_parsed_csv['temperature']
            target_press = k_target_values_parsed_csv['pressure']
            reactions_in_cti_file = gas.reaction_equations()
            for i,reaction in enumerate(target_reactions): 
                #ask about the mixture composition
                if target_press[i] == 0:
                    pressure = 1e-9
                else:
                    pressure = target_press[i]
                    
                gas.TPX = target_temp[i],pressure*101325,{'Ar':.99}
                reaction_number_in_cti = reactions_in_cti_file.index(reaction)
                k = gas.forward_rate_constants[reaction_number_in_cti]
                idx = reaction_list_from_mechanism.index(reaction)
                target_value_ks[idx].append(k)
                #check and make sure we are subtracting in the correct order 

            return target_value_ks
    
    
        if bool(target_value_rate_constant_csv):
            
            #make two unique
            unique_reactions_optimized=[]
            unique_reactions_original = []
            
            reaction_list_from_mechanism_original = gas_original.reaction_equations()
            reaction_list_from_mechanism = gas_optimized.reaction_equations()
            k_target_value_csv = pd.read_csv(target_value_rate_constant_csv)     
            for row in range(k_target_value_csv.shape[0]):
                unique_reactions_optimized.append(reaction_list_from_mechanism.index(k_target_value_csv['Reaction'][row]))
                unique_reactions_original.append(reaction_list_from_mechanism_original.index(k_target_value_csv['Reaction'][row]))
            unique_reactions_optimized = unique_list(unique_reactions_optimized)
            unique_reactions_original = unique_list(unique_reactions_original)
            
            sigma_list_for_target_ks_optimized = (self.S_matrix,k_target_value_csv,unique_reactions_optimized,gas_optimized,self.covarience)
            sigma_list_for_target_ks_original = (self.S_matrix_original,k_target_value_csv,unique_reactions_original,gas_original,self.original_covariance)
          ######################  
            target_value_temps_optimized,target_value_ks_optimized = sort_rate_constant_target_values(k_target_value_csv,unique_reactions_optimized,gas_optimized)
            target_value_temps_original,target_value_ks_original = sort_rate_constant_target_values(k_target_value_csv,unique_reactions_original,gas_original)
           ############################################# 
            target_value_ks_calculated_with_cantera_optimized = calculating_target_value_ks_from_cantera_for_sigmas(k_target_value_csv,gas_optimized,unique_reactions_optimized)
            target_value_ks_calculated_with_cantera_original = calculating_target_value_ks_from_cantera_for_sigmas(k_target_value_csv,gas_original,unique_reactions_original)
           
            for i,reaction in enumerate(unique_reactions_optimized):
                plt.figure()
                Temp_optimized,k_optimized = rate_constant_over_temperature_range_from_cantera(reaction,
                                                                  gas_original,
                                                                  initial_temperature=250,
                                                                  final_temperature=2500,
                                                                  pressure=1,
                                                                  conditions={'H2':2,'O2':1,'Ar':4})
                
                plt.semilogy(Temp_optimized,k_optimized,'b')
                #calculate sigmas 
                high_error_optimized = np.exp(sigma_list_for_target_ks_optimized[i])
                high_error_optimized = np.multiply(high_error_optimized,target_value_ks_calculated_with_cantera_optimized[i])
                low_error_optimized = np.exp(sigma_list_for_target_ks_optimized[i]*-1)
                low_error_optimized = np.multiply(low_error_optimized,target_value_ks_calculated_with_cantera_optimized[i])    
                
                plt.semilogy(target_value_temps_optimized[i],high_error_optimized,'b','--')    
                plt.semilogy(target_value_temps_optimized[i],low_error_optimized,'b','--')
                    
                    
                Temp_original,k_original = rate_constant_over_temperature_range_from_cantera(unique_reactions_original[unique_reactions_original.index(reaction)],
                                                                  gas_original,
                                                                  initial_temperature=250,
                                                                  final_temperature=2500,
                                                                  pressure=1,
                                                                  conditions={'H2':2,'O2':1,'Ar':4})
                
                plt.semilogy(Temp_original,k_original,'r')
                
                high_error_original = np.exp(sigma_list_for_target_ks_original[unique_reactions_original.index(reaction)])
                high_error_original = np.multiply(high_error_original,target_value_ks_calculated_with_cantera_original[unique_reactions_original.index(reaction)])
                
                
                low_error_original = np.exp(sigma_list_for_target_ks_original[unique_reactions_original.index(reaction)]*-1)
                low_error_original = np.multiply(low_error_original,target_value_ks_calculated_with_cantera_original[unique_reactions_original.index(reaction)])  
                
                plt.semilogy(target_value_temps_original[unique_reactions_original.index(reaction)],high_error_original,'r','--')
                plt.semilogy(target_value_temps_original[unique_reactions_original.index(reaction)],low_error_original,'r','--')
                plt.semilogy(target_value_temps_optimized[i],target_value_ks_optimized[i],'o',color='black')
                
                #add lables and titles to graphs 
                
                
                
        



            
            

                
                

                
                
            
        
                    
                    
                    
                    
                    
                    
                    
                    
