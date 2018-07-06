from ..cti_core import cti_processor as ctp

class Simulation(object):
    pasc_to_atm = 101325
    def __init__(self,pressure:float,temperature:float,observables:list,kineticSens:int,physicalSens:int
            ,conditions:dict,processor:ctp.Processor=None,cti_path = ""):
        '''
        Input:
            - pressure = float, pressure in [atm]
            - temperature = float, temperature in [K]
            - observables = list, species which sensitivity analysis is perfomred for
            - kineticSen = integer, 0 for off, 1 for on 
            - physicalSens = integer, 0 for off, 1 for on 
            - processor = ctp.Processor
            - cti_path = string, path to cti file, will construct an internal processor
            
        '''
        if processor!=None and cti_path!="":
            print("Error: Cannot give both a processor and a cti file path, pick one")
        elif processor==None and cti_path=="":
            print("Error: Must give either a processor or a cti file path")
        if processor != None:
            self.processor = processor 
        elif cti_path!="":
            self.processor = ctp.Processor(cti_path)
        self.pressure = pressure
        self.temperature = temperature
        self.observables = observables
        self.kineticSens = kineticSens
        self.physicalSens = physicalSens
        self.conditions = conditions
        
    def setTPX(self,temperature:float=-1,pressure:float=-1,conditions:dict={}):
        '''
        Set solution object for a simulation
        '''
        if temperature== -1:
            temperature = self.temperature
        if pressure == -1:
            pressure = self.pressure
        if conditions == {}:
            conditions = self.conditions
        else:
            conditions_copy = self.conditions
            for x in conditions.keys():
                conditions_copy[x] = conditions_copy[x]+conditions[x]
            conditions = conditions_copy 
        self.processor.solution.TPX=temperature,pressure*self.pasc_to_atm,conditions
    
    #always overwritten since each simulation is very different
    def run(self):
        print("Error: Simulation class itself does not implement the run method, please run a child class")


    def sensitivity_adjustment(self,temp_del:float=0.0, pres_del:float=0.0, spec_del:(str,float)=('',0.0)):
        self.setTPX(self.temperature+temp_del,
                    self.pressure+pres_del,
                    {spec_del[0]:spec_del[1]})
        #self.speciesSensitivty+=spec_del, not this easy
        data = self.run()
        self.setTPX()
        return data
    

    def species_adjustment(self,spec_del:float=0.0):
    inert_species=['Ar','AR','HE','He','Kr','KR',
                   'Xe','XE','NE','Ne']
        for x in self.conditions.keys():
            if x not in inert_species:
                sensitivity_adjustment(spec_del=(x,spec_del))

'''    #integrate with sens adjustment 
    def species_adjustment(spec_del={}):
        gas = self.solutionObject()
        ar = []
        for x in np.nditer(gas.TPX[2]):
            if x != 0:
                temp1 = np.where(gas.TPX[2] == x)
                temp2 = temp1[0]
                if np.shape(temp2)[0] == 1:
                    ar.append(temp2[0])
                    
                else:
                    for x in np.nditer(temp2):
                        if x not in ar:
                            ar.append(int(x))
        temp3 = gas.TPX[2]
        timeHistoryList = []
        for x in ar:
            temp3 = gas.TPX[2]
            temp3[x] = self.parameterAdjustment(temp3[x])
            gas.TPX = self.temperature,self.pressure*101325,self.conditions
            timeHistory = self.shockTubeSimulation()
            timeHistoryList.append(timeHistory)
            
            
       

        return timeHistoryList '''
