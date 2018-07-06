import itertools
import numpy as np
import cantera as ct
import pandas as pd
import re
from .. import simulation as sim
from ...cti_core import cti_processor as ctp

class shockTube(sim.Simulation):
    
    def __init__(self,pressure:float,temperature:float,observables:list,kineticSens:int,physicalSens:int
            ,conditions:dict,initialTime,finalTime,thermalBoundary,mechanicalBoundary,processor:ctp.Processor=None,
            cti_path=""):

        '''
        Child class of shock Tubes. Inherits all attributes and
        methods including __init__(). Also has its own internal
        init method due to additional data requirements
    
        Input:
            - initialTime = float, time the simulation will begin at in [s]
            - finalTime = float, time the simulation will end at in [s]
            - thermalBoundary = string, the boundary condtion for the shocktube.
              For example, adiabatic, or isothermal
            - mechanicalBoundary = string, the thermal boundary condition for
              the shocktube. For example, constant pressure or constant volume
        '''
        sim.Simulation.__init__(self,pressure,temperature,observables,kineticSens,physicalSens,
                                conditions,processor,cti_path)
        self.initialTime = initialTime
        self.finalTime = finalTime
        self.thermalBoundary = thermalBoundary
        self.mechanicalBoundary = mechanicalBoundary
        self.kineticSensitivities= None
        self.timeHistory = None
        self.setTPX()
    def printVars(self):
        print('initial time: {0}\nfinal time: {1}\n'.format(self.initialTime,self.finalTime),
              '\nthermalBoundary: {0}\nmechanicalBoundary: {1}'.format(self.thermalBoundary,self.mechanicalBoundary),
              '\npressure: {0}\ntemperature: {1}\nobservables: {2}'.format(self.pressure,self.temperature,self.observables),
              '\nkineticSens: {0}\nphysicalSens: {1}'.format(self.kineticSens,self.physicalSens),
              '\nTPX: {0}'.format(self.processor.solution.TPX)
              )
    def settingShockTubeConditions(self):
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

    def sensitivity_adjustment(self,temp_del:float=0.0,
                               pres_del:float=0.0,
                               spec_pair:(str,float)=('',0.0)):
        kin_temp = self.kineticSens
        self.kineticSens = 0
        data = sim.Simulation.sensitivity_adjustment(self,temp_del,pres_del,spec_pair)
        self.kineticSens = kin_temp
        return data

    def run(self,initialTime:float=-1.0, finalTime:float=-1.0):
        if initialTime == -1.0:
            initialTime = self.initialTime 
        if finalTime == -1.0:
            finalTime = self.finalTime
        self.timeHistory = None
        self.kineticSensitivities= None
        conditions = self.settingShockTubeConditions()
        mechanicalBoundary = conditions[1]
        
        #same solution for both cp and cv sims
        if mechanicalBoundary == 'constant pressure':
            shockTube = ct.IdealGasConstPressureReactor(self.processor.solution,
                                                        name = 'R1',
                                                        energy = conditions[0])
        else:
            shockTube = ct.IdealGasReactor(self.processor.solution,
                                           name = 'R1',
                                           energy = conditions[0])
        sim = ct.ReactorNet([shockTube])

        columnNames = [shockTube.component_name(item) for item in range(shockTube.n_vars)]
        columnNames = ['time']+['pressure']+columnNames
        self.timeHistory = pd.DataFrame(columns=columnNames)

        if self.kineticSens == 1:
            for i in range(self.processor.solution.n_reactions):
                shockTube.add_sensitivity_reaction(i)
            dfs = [pd.DataFrame() for x in range(len(self.observables))]
            tempArray = [np.zeros(self.processor.solution.n_reactions) for x in range(len(self.observables))]

        t = self.initialTime
        counter = 0
        while t < self.finalTime:
            t = sim.step()
            state = np.hstack([t,shockTube.thermo.P, shockTube.mass,
                               shockTube.T, shockTube.thermo.X])
            self.timeHistory.loc[counter] = state
            if self.kineticSens == 1:
                counter_1 = 0
                for observable,reaction in itertools.product(self.observables, range(self.processor.solution.n_reactions)):
                    tempArray[self.observables.index(observable)][reaction] = sim.sensitivity(observable,
                                                                                                    reaction)
                    counter_1 +=1
                    if counter_1 % self.processor.solution.n_reactions == 0:
                        dfs[self.observables.index(observable)] = dfs[self.observables.index(observable)].append(((
                            pd.DataFrame(tempArray[self.observables.index(observable)])).transpose()),
                            ignore_index=True)
            counter+=1

        if self.kineticSens == 1:
            numpyMatrixsksens = [dfs[dataframe].values for dataframe in range(len(dfs))]
            self.kineticSensitivities = np.dstack(numpyMatrixsksens)
            return self.timeHistory,self.kineticSensitivities
        else:
            return self.timeHistory
