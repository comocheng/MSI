import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import MSI.simulations.absorbance.curve_superimpose as csp 
import MSI.simulations.yaml_parser as yp
import cantera as ct

yaml_file_list = [('MSI/data/test_data/Hong_4.yaml','MSI/data/test_data/Hong_4_abs.yaml'),('MSI/data/test_data/Hong_4.yaml',)]
yaml_instance = yp.Parser()
list_of_yaml_objects = yaml_instance.load_yaml_list(yaml_list=yaml_file_list)
list_of_experiment_dicts = yaml_instance.parsing_multiple_dictonaries(list_of_yaml_objects = list_of_yaml_objects)

