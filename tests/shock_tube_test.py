import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import cantera as ct
import matplotlib.pyplot as plt

test_p = pr.Processor('MSI/data/test_data/hong_reproduction_of_data.cti')
test_tube = st.shockTube(pressure=1.653,
                         temperature=1283,
                         observables=['OH','H2O'],
                         kineticSens=0,
                         physicalSens=0,
                         conditions={'H2O2':0.003094,'H2O':0.001113,'O2':0.000556,'Ar':0.995237},
                         initialTime=0,
                         finalTime=0.001,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant pressure',
                         processor=test_p,
                         save_timeHistories=1,
                         save_physSensHistories=1)
test_tube.run()
test_tube.printVars()
time_History = test_tube.timeHistory


plt.plot(time_History['time']*1e3,time_History['OH']*1e6)
time_OH = time_History['time']
OH_ppm = time_History['OH']*1e6
plt.figure()
plt.plot(time_History['time']*1e3,time_History['H2O']*1e6)
time_H2O = time_History['time']
H2O_ppm = time_History['H2O']*1e6
