import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import MSI.simulations.absorbance.curve_superimpose as csp 
import MSI.simulations.yaml_parser as yp
import cantera as ct


test_p = pr.Processor('MSI/data/test_data/chem.cti')
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

test_tube.run()

parser = yp.Parser()
exp_loaded = parser.load_to_obj('MSI/data/test_data/Troe_6.yaml')
abs_loaded = parser.load_to_obj('MSI/data/test_data/Troe_6_abs.yaml')
abs_data = csp.superimpose_shock_tube([],test_tube.timeHistory,abs_loaded,30,[])
print(abs_data)
#loaded_tube = parser.parse_shock_tube_obj(loaded_exp=exp, loaded_absorption=absp)
#uneeded for just testing absorbance

