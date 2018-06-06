import itertools
import numpy as np
import cantera as ct
import pandas as pd
import ..simulations.simulations as sim

class ShockTube(object):
    atm_pasc = 101325
    
    def __init__(self, ctiFile, speciesNames, pressure, temperature,
            conditions, initialTime, finalTime, thermalBoundary, observalbles[],
            physical_params=[]
