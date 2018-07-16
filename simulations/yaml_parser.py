import yaml

# subpackage for reading yaml files that describe simulations and absorbance data
class Parser(object):
    def __init__(self, exp_file_path:str = '', absorption_file_path:str =  ''):
        '''
        Input:
           exp_file_path: path to yaml file which contains experimental information
           absorbtion_file_path: path to absorbtion yaml file which contains absorbance data
        '''
        self.exp_file_path = exp_file_path 
        self.absorbtion_file_path = absorption_file_path

    #config is a dict containing the yaml information
    def load_to_obj(self, path:str = ''):
        with open(path) as f:
            config = yaml.load(f)
        return config
    
    def parse_shock_tube_obj(self,loaded_exp:dict={}, loaded_absorbtion:dict={}):
     
        pressure = loaded_exp['common-properties']['pressure']['value']
        temperature = loaded_exp['common-properties']['temperature']['value']
        mole_fractions = [((concentration['mole-fraction'])) for concentration in loaded_exp['common-properties']['composition']]
        mole_fractions = [float(elm) for elm in molef_fractions]
        species_names = [(species['species']) for species in loaded_exp['common-properties']['composition']]
        conditions = dict(zip(speciesNames,moleFractions))
        thermal_boundary = loaded_exp['common-properties']['assumptions']['thermal-boundary']
        mole_fraction_observables = [point['targets'][0]['name'] for point in loaded_exp['datapoints']['mole-fraction']]
        species_uncertainties = [uncert['relative-uncertainty'] for uncert in loaded_exp['common-properties']['composition']]
        
        concentration_observables = [datapoint['targets'][0]['name'] for datapoint in loaded_exp['datapoints']['concentration']]            
        observables = [x for x in (mole_fraction_observables + concentration_observables) if x is not None]
        
        initial_time = loaded_exp['common-properties']['time']['initial-time']['value']
        #eventually going to get this from a csv file 
        final_time = loaded_exp['common-properties']['time']['final-time']['value']
    
   
        mole_fraction_csv_files = [csvfile['csvfile'] for csvfile in loaded_exp['datapoints']['mole-fraction']]
        concentration_csv_files = [csvfile['csvfile'] for csvfile in loaded_exp['datapoints']['concentration']]
        path_length = loaded_exp['apparatus']['inner-diameter']['value']
        csv_files = [x for x in (mole_fraction_csv_files + concentration_csv_files) if x is not None]


        #importing unceratinty values 
        temp_relative_uncertainty = loaded_exp['common-properties']['temperature']['relative-uncertainty']
        pressure_relative_uncertainty = loaded_exp['common-properties']['pressure']['relative-uncertainty']
        time_shift_uncertainty = loaded_exp['common-properties']['time-shift']['absolute-uncertainty']['value']
        concentration_absolute_uncertainty = [point['targets'][0]['absolute-uncertainty'] for point in loaded_exp['datapoints']['concentration']]
        concentration_relative_uncertainity = [point['targets'][0]['relative-uncertainty'] for point in loaded_exp['datapoints']['concentration']]
        mole_fraction_absolute_uncertainty = [point['targets'][0]['absolute-uncertainty'] for point in loaded_exp['datapoints']['mole-fraction']]
        mole_fraction_relative_uncertainty = [point['targets'][0]['relative-uncertainty'] for point in loaded_exp['datapoints']['mole-fraction']]        


        if loaded_absorption == {}:
            return{
               'pressure':pressure,
               'temperature':temperature,
               'conditions':conditions,
               'speciesUncertaintys':species_uncertainties,
               'thermalBoundary':thermal_boundary,
               'moleFractionObservables':mole_fraction_observables,
               'concentrationObservables': concentration_observables,               
               'observables':observables,
               'initialTime':initial_time,
               'finalTime':final_time,
               'speciesNames':species_names,
               'pathLength':path_length,
               'MoleFractions':mole_fractions,
               'moleFractionCsvFiles':mole_fraction_csv_files,
               'concentrationCsvFiles':concentration_csv_files,
               'tempRelativeUncertainty':temp_relative_uncertainty,
               'pressureRelativeUncertainty': pressure_relative_uncertainty,
               'timeShiftUncertainty':time_shift_uncertainty,
               'concentrationAbsoluteUncertainty':concentration_absolute_uncertainty,
               'concentrationRelativeUncertainity':concentration_relative_uncertainity,
               'moleFractionAbsoluteUncertainty':mole_fraction_absolute_uncertainty,
               'moleFractionRelativeUncertainty':mole_fraction_relative_uncertainty,
               'csvFiles': csv_files
           }
        
        else: #absorbtion file given
            absorbance_absolute_uncertainty = [point['absolute-uncertainty'] for point in loaded_exp['datapoints']['absorbance']]
            absorbance_relative_uncertainty = [point['relative-uncertainty'] for point in loaded_exp['datapoints']['absorbance']]
            #importing absorbance uncertainty 

            absorbance_csv_files = [csvfile['csvfile'] for csvfile in loaded_exp['datapoints']['absorbance']]
            absorbance_csv_wavelengths = [csvfile['wavelength']['value'] for csvfile in loaded_exp['datapoints']['absorbance']]
            absorption_observables = [species['species'] for species in loaded_absorption['Absorption-coefficients']]

            observables = [x for x in (mole_fraction_observables + concentration_observables + absorption_observables) if x is not None]


            uncertainty_parameter_ones = [[] for i in range(len(loaded_absorbtion['Absorption-coefficients']))]
            for uncertainty in xrange(len(loaded_absorbtion['Absorption-coefficients'])):
                temp = [wavelength['parameter-one']['absolute-uncertainty']['value'] for wavelength in loaded_absorbtion['Absorption-coefficients'][uncertainty]['wave-lengths']]
                uncertainty_parameter_ones[uncertainty] = temp
                
            uncertainty_parameter_twos = [[] for i in range(len(loaded_absorbtion['Absorption-coefficients']))]
            for uncertainty in range(len(loaded_absorbtion['Absorption-coefficients'])):
                temp = [wavelength['parameter-two']['absolute-uncertainty']['value'] for wavelength in loaded_absorption['Absorption-coefficients'][uncertainty]['wave-lengths']]
                uncertainty_parameter_twos[uncertainty] = temp            
 
            return {
                   'pressure':pressure,
                   'temperature':temperature,
                   'conditions':conditions,
                   'thermalBoundary':thermalBoundary,
                   'speciesNames': speciesNames,
                   'observables': observables,
                   'moleFractionObservables':mole_fraction_observables,
                   'concentrationObservables':concentration_observables, 
                   'absorbanceObservables':absorption_observables,
                   'initialTime': initialTime,
                   'finalTime':finalTime,
                   'speciesNames': speciesNames,
                   'MoleFractions':moleFractions,
                   'absorbanceCsvFiles': absorbance_csv_files,
                   'moleFractionCsvFiles':mole_fraction_csv_files,
                   'concentrationCsvFiles':concentration_csv_files,
                   'absorbanceCsvWavelengths': absorbance_csv_wavelengths,
                   'pathLength':path_length,
                   'tempRelativeUncertainty': temp_relative_uncertainty,
                   'pressureRelativeUncertainty': pressure_relative_uncertainty,
                   'speciesUncertaintys': species_uncertainties,
                   'timeShiftUncertainty': time_shift_uncertainty,
                   'concentrationAbsoluteUncertainty': concentration_absolute_uncertainty,
                   'concentrationRelativeUncertainity': concentration_relative_uncertainity,
                   'moleFractionAbsoluteUncertainty': mole_fraction_absolute_uncertainty,
                   'moleFractionRelativeUncertainty': mole_fraction_relative_uncertainty,
                   'absorbanceAbsoluteUncertainty': absorbance_absolute_uncertainty,
                   'absorbanceRelativeUncertainty': absorbance_relative_uncertainty,
                   'uncertaintyParameterOnes':uncertainty_parameter_ones,
                   'uncertaintyParameterTwos':uncertainty_parameter_twos,
                   }
