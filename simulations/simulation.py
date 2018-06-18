
class simulations(object):
    pasc_to_atm = 101325
    def __init__(self,ctiFile,pressure,temperature,observables,kineticSens,physicalSens):
        '''
        Input:
            - pressure = float, pressure in [atm]
            - temperature = float, temperature in [K]
            - observables = list, species which sensitivity analysis is perfomred for
            - kineticSen = integer, 0 for off, 1 for on 
            - physicalSens = integer, 0 for off, 1 for on 
            
        '''
        self.ctiFile = ctiFile
        self.pressure = pressure
        self.temperature = temperature
        self.observables = observables
        self.kineticSens = kineticSens
        self.physicalSens = physicalSens
        
    def solutionObject(self):
        '''
        Set solution object for a simulation
        '''
        gas = ct.Solution(self.ctiFile)
        gas.TPX = self.temperature,self.pressure*pasc_to_atm,self.conditions
        
        return gas
    
    #always overwritten since each simulation is very different
    #possible do have a main loop, but call run in this loop?
    def run(self, initialTime, finalTime):
        return (initialTime,finalTime)
