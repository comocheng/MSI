import MSI.simulations as sim
import re

def superimpose_shock_tube(absorbance_csv_files:list,
                           simulation:sim.instruments.shock_tube.shockTube,
                           absorb:dict,pathlength:float,
                           absorbance_csv_wavelengths:list):
    '''Input:
        absorbanceCsvFiles: list of absorbance experimental data
        simulation: shockTube instance for accessing TPX and solution data
        absorb: dict of yaml parsed absorbance data ie yaml_parser.parse_shock_tube_obj was run
        pathLenth: diameter of shock tube, integer
        absorbanceCsvWavelengths: list of wavelengths absorbances were measured at, in nm
    '''
    print('Importing shock tube absorbance data the following csv files...') 
    print(absorbance_csv_files)
    print(absorbance_csv_wavelengths)
    
    coupled_coefficients = couple_parameters(absorb)
    if coupled_coefficients == -1:
        print("Error: could not construct coupled coefficients")
        return -1

    wavelengths = get_wavelengths(absorb)
    #functional form takes A B C D, changes which equation used
    functional_form = get_funtional(absorb) 
    
    #group data by species for easier manipulation
    species_and_wavelengths = dict(zip(species, wavelengths))
    species_and_coupled_coefficients = dict(zip(species,coupled_coefficients))
    species_and_functional_form = dict(zip(species,functional_form)) 
    
    temperature_matrix = simulation.processor.solution['temperature'].as_matrix()

    flat_list = [item for sublist in wavelengths for item in sublist]
    flat_list = list(set(flat_list))
    absorbance_species_list = []
    for i, wl in enumerate(flat_list):
        absorbance_species.append([])
        for j,specie in enumerate(species):
            if wl in SandWL[specie]:
                absorbance_species_list[i].append(specie)
 
    d = {}
    s = {}
    ATemp= {}
    APressure = {}
    ASpecies = {}

    absorbance_species_wavelengths= []
    for i in xrange(len(absorbanceSpeciesList)):
        wavelength = flat_list[i]
        for j in xrange(len(absorbanceSpeciesList[i])):            
            value = absorbanceSpeciesList[i][j]
            index = species_and_wavelengths[value].index(wavelength)
            absorbance_species_wavelengths.append(calc_absorb(simulation,
                                                          value,
                                                          species_and_functional_form[value][index],
                                                          species_and_coupled_coefficients[value][index],
                                                          wavelength,
                                                          pathlength),
                                              value,
                                              wavelength)

    return absorbance_species_wavelengths

def calc_absorb(simulation:sim.instruments.shock_tube.shockTube,
                   species,
                   ff,
                   cc,
                   wavelength:float,
                   pathlength):
    if ff == 'A':
        epsilon = ((cc[1]*Temperature) + cc[0])
    if ff == 'B':
        epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))*1000
    if ff == 'C':
        epsilon = cc[0] 
    
    if wavelength == 215: #does this really need to be here?
       epsilon *= 1000
       #multiplying by 1000 to convert from L to cm^3 from the epsilon given in paper 
       #this applies if the units on epsilon are given as they are in kappl paper 
       #must calcuate and pass in reactor volume 
    concentration = ((np.true_divide(1,simulation.processor.solution['temperature'].as_matrix().flatten()))*(simulation.processor.solution['pressure'].as_matrix().flatten()) * (1/(8.314e6)))*simulation.processor.solution[value].as_matrix().flatten()
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
        temp = [wlh['parameter-one']['value'] for wl in absorb['Absorption-coefficients'][p1]['wave-lengths']]
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

    coupled_coefficients = [zip(parameter_ones[x],parameter_twos[x]) for x in range(len(parameter_ones))]
    return coupled_coefficients 


def get_funtional(absorb:dict):
    functional_form = []
    for form in range(len(absorb['Absorption-coefficients'])):
        temp = [wl['functional-form'] for wl in absorb['Absorption-coefficients'][form]['wave-lengths']]
        functional_form.append(temp)
    return functional_form
