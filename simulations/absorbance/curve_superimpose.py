import MSI.simulations as sim
import numpy as np

class Absorb:
    def __init__(self):
        self.saved_abs_data = []
        self.saved_perturb_data = [] 
    def perturb_abs(self,del_param:float,simulation:sim.instruments.shock_tube.shockTube,
                    absorb:dict,pathlength:float,
                    absorbance_csv_files:list=[],
                    absorbance_csv_wavelengths:list=[],
                    kinetic_sens=0):

        #modify absorb dict, modify each coef
        #remember to put it back
        #essentially modify the coupled coefficient before calculating the absorbance data
        coupled_coefficients = self.couple_parameters(absorb) 
        #loop over the coupled_coefficients, replace, then put back
        self.saved_perturb_data.clear()
        species = [species['species'] for species in absorb['Absorption-coefficients']]
        species_and_coupled_coefficients = dict(list(zip(species,coupled_coefficients)))
        species_and_wavelengths = dict(list(zip(species, self.get_wavelengths(absorb))))
        for species in species_and_coupled_coefficients.keys():
            index = 0
            for cc in species_and_coupled_coefficients[species]:
                orig_cc = cc
                
                cc = (cc[0]+del_param,cc[1])
                changed_data = self.get_abs_data(simulation,
                                                 absorb,
                                                 pathlength,
                                                 kinetic_sens = 1)
                self.saved_perturb_data.append(((cc,species,species_and_wavelengths[species][index]),
                                                changed_data))
                
                cc = (orig_cc[0],cc[1]+del_param)
                changed_data = self.get_abs_data(simulation,
                                                 absorb,
                                                 pathlength,
                                                 kinetic_sens = 1)
                self.saved_perturb_data.append(((species,species_and_wavelengths[species][index],cc),
                                                changed_data))
                cc = orig_cc
                index += 1
        return self.saved_perturb_data
    
    def superimpose_shock_tube(self,simulation:sim.instruments.shock_tube.shockTube,
                               absorb:dict,pathlength:float,
                               absorbance_csv_files:list=[],
                               absorbance_csv_wavelengths:list=[],
                               kinetic_sens=0):
        '''Input:
            absorbanceCsvFiles: list of absorbance experimental data
            time_history: time history from some run used to get temperature matrix
            absorb: dict of yaml parsed absorbance data ie yaml_parser.parse_shock_tube_obj was run
            pathLenth: diameter of shock tube, integer
            absorbanceCsvWavelengths: list of wavelengths absorbances were measured at, in nm
        '''
        print('Importing shock tube absorbance data the following csv files...') 
        #print(absorbance_csv_files)
        #print(absorbance_csv_wavelengths)
        
        
        abs_data = self.get_abs_data(simulation,
                                absorb,
                                pathlength,
                                kinetic_sens)
        self.saved_abs_data.append(abs_data)
        return abs_data 
    
    def get_abs_data(self,simulation,absorb,pathlength,kinetic_sens = 0):
        
        coupled_coefficients = self.couple_parameters(absorb)
        if coupled_coefficients == -1:
            print("Error: could not construct coupled coefficients")
            return -1
        
        species = [species['species'] for species in absorb['Absorption-coefficients']]
        wavelengths = self.get_wavelengths(absorb)
        #functional form takes A B C D, changes which equation used
        functional_form = self.get_functional(absorb) 
        #group data by species for easier manipulation
        species_and_wavelengths = dict(list(zip(species, wavelengths)))
        species_and_coupled_coefficients = dict(list(zip(species,coupled_coefficients)))
        species_and_functional_form = dict(list(zip(species,functional_form))) 
        
        flat_list = [item for sublist in wavelengths for item in sublist]
        flat_list = list(set(flat_list))
        absorbance_species_list = []
        for i, wl in enumerate(flat_list):
            absorbance_species_list.append([])
            for j,specie in enumerate(species):
                if wl in species_and_wavelengths[specie]:
                    absorbance_species_list[i].append(specie)
     
        absorbance_species_wavelengths= []
        for i in range(len(absorbance_species_list)):
            wavelength = flat_list[i]
            for j in range(len(absorbance_species_list[i])):            
                value = absorbance_species_list[i][j]
                index = species_and_wavelengths[value].index(wavelength)
                absorbance_species_wavelengths.append((self.calc_absorb(value,
                                                                  species_and_functional_form[value][index],
                                                                  species_and_coupled_coefficients[value][index],
                                                                  wavelength,
                                                                  pathlength,
                                                                  simulation.timeHistories[0]),
                                                       value,
                                                       wavelength))

        
        summed_data = {} 
        for x in absorbance_species_wavelengths:
            if x[2] not in summed_data.keys():
                summed_data[x[2]] = x[0]
            else:
                summed_data[x[2]] += x[0]
            
        if kinetic_sens == 0:
            return summed_data
        else:
            return summed_data, self.calc_abs_sens(simulation, 
                                              species_and_wavelengths,
                                              species_and_functional_form,
                                              species_and_coupled_coefficients,
                                              absorbance_species_wavelengths,
                                              pathlength,
                                              summed_data)

    def calc_abs_sens(self,simulation,
                      species_and_wavelengths,
                      species_and_functional_form,
                      species_and_coupled_coefficients,
                      absorbance_species_wavelengths,
                      pathlength,
                      summed_absorbtion):
        species_and_sensitivities = {}
        for x,i in enumerate(simulation.observables):
            slice_2d = simulation.kineticSensitivities[:,:,x]
            species_and_sensitivities[i]=slice_2d
        temperature_matrix = simulation.timeHistories[0]['temperature'].values
        pressure_matrix = simulation.timeHistories[0]['pressure'].values

        ind_wl_derivs = {} 
        for species in species_and_sensitivities.keys():
            if species not in species_and_wavelengths.keys():
                continue
            #do epsilon and con calc, then mult
            wavelengths = species_and_wavelengths[species]
            
            for j in range(0,len(wavelengths)):
                wavelength = species_and_wavelengths[species][j] 
                if wavelength not in ind_wl_derivs.keys():
                    net_sum = np.zeros(shape=(simulation.kineticSensitivities.shape[0:2])) #only need 2d info, since sum over observables
                else:
                    net_sum = ind_wl_derivs[wavelength]
                index = species_and_wavelengths[species].index(wavelength)
                cc = species_and_coupled_coefficients[species][index]
                ff = species_and_functional_form[species][index]
                if ff == 'A':
                    epsilon = ((cc[1]*temperature_matrix) + cc[0])
                if ff == 'B':
                    epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],temperature_matrix)))))
                if ff == 'C':
                    epsilon = cc[0] 
                
                if wavelength == 215: #does this really need to be here, takes care of specific paper case?
                   epsilon *= 1000
                
                concentration = np.true_divide(1,temperature_matrix.flatten())*pressure_matrix.flatten()
                concentration *= (1/(8.314e6))*simulation.timeHistories[0][species].values.flatten()
                temp = np.multiply(species_and_sensitivities[species],concentration.reshape((np.shape(concentration)[0],1)))
                net_sum += np.multiply(temp,epsilon.reshape((np.shape(epsilon)[0],1)))
                
                ind_wl_derivs[wavelength]=net_sum
        
        for single_wl in ind_wl_derivs.keys():
            flat_list = np.array(list(summed_absorbtion[single_wl]))
            for i in range(0,len(flat_list)):
                flat_list[i] = 1/flat_list[i]
            for column in ind_wl_derivs[single_wl].T:
                column*=flat_list.flatten()
        
        return ind_wl_derivs


    def calc_absorb(self,species,
                    ff,
                    cc,
                    wavelength,
                    pathlength,
                    time_history):
        
        temperature_matrix = time_history['temperature'].values
        pressure_matrix = time_history['pressure'].values
        if ff == 'A':
            epsilon = ((cc[1]*temperature_matrix) + cc[0])
        if ff == 'B':
            epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],temperature_matrix)))))
        if ff == 'C':
            epsilon = cc[0] 
        
        if wavelength == 215: #does this really need to be here?
           epsilon *= 1000
           #multiplying by 1000 to convert from L to cm^3 from the epsilon given in paper 
           #this applies if the units on epsilon are given as they are in kappl paper 
           #must calcuate and pass in reactor volume 
        concentration = np.true_divide(1,temperature_matrix.flatten())*pressure_matrix.flatten()
        concentration *= (1/(8.314e6))*time_history[species].values.flatten()
        
        
        absorb = pathlength*(epsilon*concentration)
        return absorb
     
    def get_wavelengths(self,absorb:dict):
        wavelengths = []  #get the wavelengths
        for sp in range(len(absorb['Absorption-coefficients'])):
            temp = [wl['value'] for wl in absorb['Absorption-coefficients'][sp]['wave-lengths']]
            wavelengths.append(temp)

        return wavelengths

    def couple_parameters(self,absorb:dict()):
        parameter_ones = [] #for the epsilon calculations, there is always a parameter
        for p1 in range(len(absorb['Absorption-coefficients'])):
            temp = [wl['parameter-one']['value'] for wl in absorb['Absorption-coefficients'][p1]['wave-lengths']]
            parameter_ones.append(temp)
        
        #parameter two does not always really exist, but we will always define it in the yaml file to be 0 in that case
        #do this for easy index matching
        parameter_twos = [] 
        for p2 in range(len(absorb['Absorption-coefficients'])):
            temp = [wl['parameter-two']['value'] for wl in absorb['Absorption-coefficients'][p2]['wave-lengths']]
            parameter_twos.append(temp)
        if len(parameter_ones) != len(parameter_twos):
            print("Error: number of parameters do not match, change the yaml file")
            return -1

        coupled_coefficients = [list(zip(parameter_ones[x],parameter_twos[x])) for x in range(len(parameter_ones))]
        return coupled_coefficients 


    def get_functional(self,absorb:dict):
        functional_form = []
        for form in range(len(absorb['Absorption-coefficients'])):
            temp = [wl['functional-form'] for wl in absorb['Absorption-coefficients'][form]['wave-lengths']]
            functional_form.append(temp)
        return functional_form
