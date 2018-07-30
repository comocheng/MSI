import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import MSI.optimization.matrix_loader as ml
import cantera as ct
import pandas
test_p = pr.Processor('MSI/data/test_data/FFCM1.cti')
test_tube = st.shockTube(pressure=1.74,
                         temperature=1880,
                         observables=['OH','H2O'],
                         kineticSens=1,
                         physicalSens=0,
                         conditions={'H2O':.013,'O2':.0099,'H':.0000007,'Ar':0.9770993},
                         initialTime=0,
                         finalTime=0.5,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant pressure',
                         processor=test_p,
                         save_timeHistories=1,
                         save_physSensHistories=1)

csv_paths = ['MSI/data/test_data/hong_h2o_4.csv','MSI/data/test_data/hong_oh_4.csv']
exp_data = test_tube.importExperimentalData(csv_paths)

test_tube.run() #set up original time history
int_ksens_exp_mapped= test_tube.map_and_interp_ksens()#ksens is wiped on rerun so int it before
test_tube.sensitivity_adjustment(temp_del = .01)
test_tube.sensitivity_adjustment(pres_del = .01)
test_tube.species_adjustment(.01) #do some sensitivity adjustments

int_tp_psen_against_experimental = test_tube.interpolate_experimental([test_tube.interpolate_physical_sensitivities(index=1),
                                                                    test_tube.interpolate_physical_sensitivities(index=2)])
int_spec_psen_against_experimental = test_tube.interpolate_experimental(pre_interpolated=test_tube.interpolate_species_sensitivities())


test_p2 = pr.Processor('MSI/data/test_data/FFCM1.cti')
test_tube2 = st.shockTube(pressure=1.74,
                         temperature=1880,
                         observables=['OH','H2O'],
                         kineticSens=1,
                         physicalSens=0,
                         conditions={'H2O':.013,'O2':.0099,'H':.0000007,'Ar':0.9770993},
                         initialTime=0,
                         finalTime=0.5,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant pressure',
                         processor=test_p,
                         save_timeHistories=1,
                         save_physSensHistories=1)

csv_paths2 = ['MSI/data/test_data/hong_h2o_4.csv','MSI/data/test_data/hong_oh_4.csv']
exp_data2 = test_tube.importExperimentalData(csv_paths)

test_tube2.run() #set up original time history
int_ksens_exp_mapped2= test_tube.map_and_interp_ksens()#ksens is wiped on rerun so int it before
test_tube.sensitivity_adjustment(temp_del = .01)
test_tube.sensitivity_adjustment(pres_del = .01)
test_tube.species_adjustment(.01) #do some sensitivity adjustments

int_tp_psen_against_experimental2 = test_tube.interpolate_experimental([test_tube2.interpolate_physical_sensitivities(index=1),test_tube2.interpolate_physical_sensitivities(index=2)])
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
