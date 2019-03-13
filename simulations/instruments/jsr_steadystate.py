# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 15:33:46 2018

@author: Mark Barbet
"""

import cantera as ct
from .. import simulation as sim
from ...cti_core import cti_processor as ctp
import pandas as pd
import numpy as np
import time

class JSR_steadystate(sim.Simulation):
    
    '''child class of sim.Simulaton.  Inherits all attributes and methods including __init__().  
    Also has internal init due to data requirements'''
    
    
    
    
    
    def __init__(self,pressure:float,temperature:float,observables:list,
                 kineticSens:int,physicalSens:int,conditions:dict,thermalBoundary,mechanicalBoundary,
                 processor:ctp.Processor=None,cti_path="", 
                 save_physSensHistories=0,moleFractionObservables:list=[],
                 absorbanceObservables:list=[],concentrationObservables:list=[],
                 fullParsedYamlFile:dict={},residence_time:float=1.0,pvalveCoefficient:float=0.01,maxpRise:float=0.001):
        
        sim.Simulation.__init__(self,pressure,temperature,observables,kineticSens,physicalSens,
                                conditions,processor,cti_path)
        self.thermalBoundary = thermalBoundary
        self.mechanicalBoundary = mechanicalBoundary
        self.kineticSensitivities= None
        self.experimentalData = None
        self.concentrationObservables = concentrationObservables
        self.moleFractionObservables = moleFractionObservables
        self.absorbanceObservables = absorbanceObservables
        self.fullParsedYamlFile =  fullParsedYamlFile
        self.pvalveCoefficient=pvalveCoefficient
        self.maxPrise=maxpRise
        self.energycon='off'
        self.residence_time=residence_time
        
        if save_physSensHistories == 1:
            self.physSensHistories = []
        self.setTPX()
        self.dk = 0.01
        self.rtol=1e-6
        self.atol=1e-6
        self.solution=None
        
        
    def set_geometry(self,volume=0.1):
        self.reactor_volume=volume
        
        
    def printVars(self):
        print()
        
    def settingJSRConditions(self):
        '''
        Determine the mechanical and thermal boundary conditions for a 
        shock tube.
        '''
        
        #assigning the thermal boundary variable
        if re.match('[aA]diabatic',self.thermalBoundary):
            energy = 'on'
        elif re.match('[iI]sothermal',self.thermalBoundary):
            energy = 'off'
        else:
            raise Exception('Please specify a thermal boundary condition, adiabatic or isothermal')
        #assigning the mehcanical boundary variable 
        if re.match('[Cc]onstant [Pp]ressure',self.mechanicalBoundary):
            mechBoundary = 'constant pressure'
        elif re.match('[Cc]onstant [Vv]olume',self.mechanicalBoundary):
            mechBoundary = 'constant volume'
        else:
            raise Exception('Please specifiy a mehcanical boundary condition, constant pressure or constant volume')
        #return the thermal and mechanical boundary of the shock tube 
        return energy,mechBoundary
    
    
    
    def run(self):
        
        gas=self.processor.solution
        reactorPressure=gas.P
        pressureValveCoefficient=self.pvalveCoefficient
        maxPressureRiseAllowed=self.maxPrise
        
        
        #Build the system components for JSR
        fuelAirMixtureTank=ct.Reservoir(gas)
        exhaust=ct.Reservoir(gas)
        
        stirredReactor=ct.IdealGasReactor(gas,energy=self.energycon,volume=self.reactor_volume)
        massFlowController=ct.MassFlowController(upstream=fuelAirMixtureTank,
                                                 downstream=stirredReactor,mdot=stirredReactor.mass/self.residence_time)
        pressureRegulator=ct.Valve(upstream=stirredReactor,downstream=exhaust,K=pressureValveCoefficient)
        reactorNetwork=ct.ReactorNet([stirredReactor])
        
        if bool(self.observables) and self.kineticSens==1:
            for i in range(gas.n_reactions):
                stirredReactor.add_sensitivity_reaction(i)
            
        if self.kineticSens and bool(self.observables)==False:
            #except:
                print('Please supply a non-empty list of observables for sensitivity analysis or set kinetic_sens=0')
        if self.physicalSens==1 and bool(self.observables)==False:
            #except:
                print('Please supply a non-empty list of observables for sensitivity analysis or set physical_sens=0')
        
        # now compile a list of all variables for which we will store data
        columnNames = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
        columnNames = ['pressure'] + columnNames

        # use the above list to create a DataFrame
        timeHistory = pd.DataFrame(columns=columnNames)

        # Start the stopwatch
        tic = time.time()
        
        
        #Establish a matrix to hold sensitivities for kinetic parameters, along with tolerances
        if self.kineticSens==1 and bool(self.observables):
            #senscolumnNames = ['Reaction']+observables     
            senscolumnNames = self.observables
            #sensArray = pd.DataFrame(columns=senscolumnNames)
            #senstempArray = np.zeros((gas.n_reactions,len(observables)))
            dfs = [pd.DataFrame() for x in range(len(self.observables))]
            tempArray = [np.zeros(self.processor.solution.n_reactions) for x in range(len(self.observables))]
            
            reactorNetwork.rtol_sensitivity = self.rtol
            reactorNetwork.atol_sensitivity = self.atol        
        
        
        reactorNetwork.advance_to_steady_state()
        final_pressure=stirredReactor.thermo.P
        sens=reactorNetwork.sensitivities()
        if self.kineticSens==1 and bool(self.observables):
            for k in range(len(self.moleFractionObservables)):
                dfs[k] = dfs[k].append(((pd.DataFrame(sens[k,:])).transpose()),ignore_index=True)
        toc = time.time()
        print('Simulation Took {:3.2f}s to compute'.format(toc-tic))
   
        columnNames = []
        #Store solution to a solution array
        #for l in np.arange(stirredReactor.n_vars):
            #columnNames.append(stirredReactor.component_name(l))
        columnNames=[stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
        #state=stirredReactor.get_state()
        state=np.hstack([stirredReactor.mass, 
                   stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X])
        data=pd.DataFrame(state).transpose()
        data.columns=columnNames
        pressureDifferential = timeHistory['pressure'].max()-timeHistory['pressure'].min()
        if(abs(pressureDifferential/reactorPressure) > maxPressureRiseAllowed):
            #except:
                print("WARNING: Non-trivial pressure rise in the reactor. Adjust K value in valve")
        
        if self.kineticSens==1:
            numpyMatrixsksens = [dfs[dataframe].values for dataframe in range(len(dfs))]
            self.kineticSensitivities = np.dstack(numpyMatrixsksens)
            self.solution=data
            return (self.solution,self.kineticSensitivities)
        else:
            return self.solution
        
        
        
        
        
class JSR_multiTemp_steadystate(sim.Simulation):
        
    def __init__(self,volume:float,pressure:float,temperatures:list,observables:list,
                 kineticSens:int,physicalSens:int,conditions:dict,thermalBoundary,mechanicalBoundary,
                 processor:ctp.Processor=None,cti_path="", 
                 save_physSensHistories=0,moleFractionObservables:list=[],
                 absorbanceObservables:list=[],concentrationObservables:list=[],
                 fullParsedYamlFile:dict={},residence_time:float=1.0,pvalveCoefficient:float=0.01,maxpRise:float=0.001):
        
#    sim.Simulation.__init__(self,pressure,temperature,observables,kineticSens,physicalSens,
#                                conditions,processor,cti_path)
        self.volume=volume
        self.temperatures=temperatures
        self.JSR_objects=[]
        for i in range(len(self.temperatures)):
            self.JSR_objects.append(JSR_steadystate(pressure,self.temperatures[i],observables,
                 kineticSens,physicalSens,conditions,thermalBoundary,mechanicalBoundary,
                 processor,cti_path, 
                 save_physSensHistories,moleFractionObservables,
                 absorbanceObservables,concentrationObservables,
                 fullParsedYamlFile,residence_time,pvalveCoefficient,maxpRise))
            
            
            
    def run(self):
        
        for i in range(len(self.JSR_objects)):
            self.JSR_objects[i].set_geometry(self.volume)
            
            self.JSR_objects[i].run()
        
        
        
        
        
        
        
        
        