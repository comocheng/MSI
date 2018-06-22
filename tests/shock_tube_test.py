import sys
sys.path.append('.')
import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import cantera as ct

test_p = pr.Processor('MSI/data/test_data/FFCM1.cti')
test_tube = st.shockTube(pressure=1.74,
                         temperature=1880,
                         observables=['OH','H2O'],
                         kineticSens=1,
                         physicalSens=0,
                         conditions={'H2O':.013,'O2':.0099,'H':.0000007,'Ar':0.9770993},
                         processor=test_p,
                         initialTime=0,
                         finalTime=0.1,
                         thermalBoundary='Adiabatic',
                         mechanicalBoundary='constant pressure')
test_tube.run()
test_tube.printVars()
print(test_tube.timeHistory)
print(test_tube.kineticSensitivities)
