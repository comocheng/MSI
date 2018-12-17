import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import cantera as ct

test_p = pr.Processor('MSI/data/test_data/FFCM1_custom_original_updated.cti')
test_tube = st.shockTube(pressure=1.74,
                         temperature=1880,
                         observables=['OH','H2O'],
                         kineticSens=1,
                         physicalSens=0,
                         conditions={'H2O':0.012438871,'O2': 0.009715602,'H':.0000007,'Ar':0.9770993},
                         initialTime=0,
                         finalTime=0.1,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant pressure',
                         processor=test_p)
test_tube.run()
test_tube.printVars()
print(test_tube.timeHistory)
print(test_tube.kineticSensitivities)
