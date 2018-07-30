import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import MSI.optimization.matrix_loader as ml
import MSI.simulations.absorbance.curve_superimpose as csp
import MSI.simulations.yaml_parser as yp
import cantera as ct
import pandas
test_p = pr.Processor('MSI/data/test_data/optimized_burke.cti')
test_tube = st.shockTube(pressure=3.44187,
                         temperature=1079,
                         observables=['H2O2','HO2','O2'],
                         kineticSens=1,
                         physicalSens=0,
                         conditions={'H2O2':0.00195373,'Ar':0.99804627},
                         initialTime=0,
                         finalTime=0.0014,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant pressure',
                         processor=test_p,
                         save_timeHistories=1,
                         save_physSensHistories=1)

csv_paths = ['MSI/data/test_data/hong_h2o_4.csv','MSI/data/test_data/hong_oh_4.csv']
exp_data = test_tube.importExperimentalData(csv_paths)

test_tube.run() #set up original time history
abs1_instance = csp.Absorb()
parser = yp.Parser()
#exp1_loaded = parser.load_to_obj('MSI/data/test_data/Troe_6.yaml')
abs1_loaded = parser.load_to_obj('MSI/data/test_data/Troe_6_abs.yaml')
abs1_data = abs1_instance.superimpose_shock_tube(test_tube,abs1_loaded,30,kinetic_sens=1)
perturbed_coef1 = abs1_instance.perturb_abs_coef(.01,
                                          test_tube,
                                          abs1_loaded,30,
                                          summed_data = abs1_data[0],
                                          dk = .01)

int_ksens_exp_mapped= test_tube.map_and_interp_ksens()#ksens is wiped on rerun so int it before
test_tube.sensitivity_adjustment(temp_del = .01)
test_tube.sensitivity_adjustment(pres_del = .01)
test_tube.species_adjustment(.01) #do some sensitivity adjustments
abs1_phys_sens = abs1_instance.absorb_phys_sensitivities(test_tube,abs1_data[0],abs1_loaded,30,dk=.01)

loaded_experimental_data1 = abs1_instance.import_experimental_data(['MSI/data/test_data/tro_6_abs_1.csv',
                                                                  'MSI/data/test_data/tro_6_abs_2.csv'])

interp_abs1_exp= abs1_instance.interpolate_experimental(test_tube,loaded_experimental_data1,
                                                        original_summed_absorption=abs1_data[0],
                                                        abs_kinetic_sens = abs1_data[1],
                                                        abs_phys_sens = abs1_phys_sens,
                                                        abs_coef_sens = perturbed_coef1)
int_tp_psen_against_experimental = test_tube.interpolate_experimental([test_tube.interpolate_physical_sensitivities(index=1),
                                                                       test_tube.interpolate_physical_sensitivities(index=2)])
int_spec_psen_against_experimental = test_tube.interpolate_experimental(pre_interpolated=test_tube.interpolate_species_sensitivities())


test_p2 = pr.Processor('MSI/data/test_data/optimized_burke.cti')
test_tube2 = st.shockTube(pressure=3.44187,
                         temperature=1079,
                         observables=['H2O2','HO2','O2'],
                         kineticSens=1,
                         physicalSens=0,
                         conditions={'H2O2':0.00195373,'Ar':0.99804627},
                         initialTime=0,
                         finalTime=0.0014,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant pressure',
                         processor=test_p,
                         save_timeHistories=1,
                         save_physSensHistories=1)

csv_paths2 = ['MSI/data/test_data/hong_h2o_4.csv','MSI/data/test_data/hong_oh_4.csv']
exp_data2 = test_tube.importExperimentalData(csv_paths)

test_tube2.run() #set up original time history
int_ksens_exp_mapped2= test_tube2.map_and_interp_ksens()#ksens is wiped on rerun so int it before
test_tube2.sensitivity_adjustment(temp_del = .01)
test_tube2.sensitivity_adjustment(pres_del = .01)
test_tube2.species_adjustment(.01) #do some sensitivity adjustments

int_tp_psen_against_experimental2 = test_tube.interpolate_experimental([test_tube2.interpolate_physical_sensitivities(index=1),
                                                                        test_tube2.interpolate_physical_sensitivities(index=2)])
int_spec_psen_against_experimental2 = test_tube2.interpolate_experimental(pre_interpolated=test_tube2.interpolate_species_sensitivities())

list_of_interpolated_kinetic_sens = [int_ksens_exp_mapped,int_ksens_exp_mapped2]
list_of_interpolated_tp_sens = [int_tp_psen_against_experimental,int_tp_psen_against_experimental2]
list_of_interpolated_species_sens = [int_spec_psen_against_experimental,int_spec_psen_against_experimental2]



mloader = ml.OptMatrix()
S = mloader.load_S(2,list_of_interpolated_kinetic_sens,
                     list_of_interpolated_tp_sens, 
                     list_of_interpolated_species_sens,
                     [],#need to do absorbance sens
                     test_tube.observables)

print(S)
