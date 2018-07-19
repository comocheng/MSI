import MSI.simulations as sim
import numpy as np
import matplotlib as mpl
def superimpose_shock_tube(absorbance_csv_files:list,
                           time_history,
                           absorb:dict,pathlength:float,
                           absorbance_csv_wavelengths:list):
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
    

    abs_data = get_abs_data(time_history,
                            absorb,
                            pathlength)
    summed_data = {} 
    for x in abs_data:
        if x[2] not in summed_data.keys():
            summed_data[x[2]] = x[0]
        else:
            summed_data[x[2]] += x[0]
    return summed_data
    
def get_abs_data(time_history,absorb,pathlength):
    
    coupled_coefficients = couple_parameters(absorb)
    if coupled_coefficients == -1:
        print("Error: could not construct coupled coefficients")
        return -1
    
    species = [species['species'] for species in absorb['Absorption-coefficients']]
    wavelengths = get_wavelengths(absorb)
    #functional form takes A B C D, changes which equation used
    functional_form = get_funtional(absorb) 
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
            absorbance_species_wavelengths.append((calc_absorb(value,
                                                              species_and_functional_form[value][index],
                                                              species_and_coupled_coefficients[value][index],
                                                              wavelength,
                                                              pathlength,
                                                              time_history),
                                              value,
                                              wavelength))

    return absorbance_species_wavelengths

def calc_absorb(species,
                   ff,
                   cc,
                   wavelength:float,
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
    concentration = ((np.true_divide(1,temperature_matrix.flatten()))*(pressure_matrix.flatten()) * (1/(8.314e6)))*time_history[species].values.flatten()
    
    absorb = pathlength*(epsilon*concentration)
    return absorb
 
def get_wavelengths(absorb:dict):
    wavelengths = []  #get the wavelengths
    for sp in range(len(absorb['Absorption-coefficients'])):
        temp = [wl['value'] for wl in absorb['Absorption-coefficients'][sp]['wave-lengths']]
        wavelengths.append(temp)

    return wavelengths

def couple_parameters(absorb:dict()):
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

    coupled_coefficients = [list(zip(parameter_ones[x],parameter_twos[x]))for x in range(len(parameter_ones))]
    return coupled_coefficients 


def get_funtional(absorb:dict):
    functional_form = []
    for form in range(len(absorb['Absorption-coefficients'])):
        temp = [wl['functional-form'] for wl in absorb['Absorption-coefficients'][form]['wave-lengths']]
        functional_form.append(temp)
    return functional_form
