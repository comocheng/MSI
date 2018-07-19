def importingAbsorbanceData(absorbanceCsvFiles, modelDataObject , absorb,  pathLength,absorbanceCsvWavelengths, speciesNames ,kinetic_sens = 0 , physical_sens=0):
    print('Importing Absorbance Data In Progress From The Following CSV Files...')    
    print(absorbanceCsvFiles)
    print(absorbanceCsvWavelengths)
    
    physical_params = ['T','P','X']
    print('these are the species names')
    print(speciesNames)
    if 'AR' in speciesNames:
        addBackin = 'AR'
        speciesNames.remove('AR')
    if 'Ar' in speciesNames:        
        addBackin = 'Ar'
        speciesNames.remove('Ar')
    if 'HE' in speciesNames:
        addBackin = 'HE'
        speciesNames.remove('HE')
    if 'He' in speciesNames:
        addBackin = 'He'
        speciesNames.remove('He')
        
    physicalParamsSpecies = ['X'+ species for species in speciesNames]
    print(physicalParamsSpecies)
    physical_params = ['T','P','X'] 
    overallList = ['T','P']+physicalParamsSpecies
    print(overallList)
    
    
    species = [species['species'] for species in absorb['Absorption-coefficients']]

    wavelengths = [[] for i in xrange(len(absorb['Absorption-coefficients']))]                      
    for specie in xrange(len(absorb['Absorption-coefficients'])):
        temp = [wavelength['value'] for wavelength in absorb['Absorption-coefficients'][specie]['wave-lengths']]
        wavelengths[specie] = temp  

    
    parameterOnes = [[] for i in xrange(len(absorb['Absorption-coefficients']))]
    for parameterOne in xrange(len(absorb['Absorption-coefficients'])):
        temp = [wavelength['parameter-one']['value'] for wavelength in absorb['Absorption-coefficients'][parameterOne]['wave-lengths']]
        parameterOnes[parameterOne] = temp
    
    parameterTwos = [[] for i in xrange(len(absorb['Absorption-coefficients']))] 
    for parameterTwo in xrange(len(absorb['Absorption-coefficients'])):
        temp = [wavelength['parameter-two']['value'] for wavelength in absorb['Absorption-coefficients'][parameterTwo]['wave-lengths']]
        parameterTwos[parameterTwo] = temp

    functionalForm = [[] for i in xrange(len(absorb['Absorption-coefficients']))]
    for form in xrange(len(absorb['Absorption-coefficients'])):
        temp = [wavelength['functional-form'] for wavelength in absorb['Absorption-coefficients'][form]['wave-lengths']]
        functionalForm[form] = temp
    
    
    SandWL = dict(zip(species, wavelengths))
    coupledCoefficients = [zip(parameterOnes[x],parameterTwos[x]) for x in xrange(len(parameterOnes))]


    
    SandCoupledCoefficients = dict(zip(species,coupledCoefficients))
    SandfunctionalForm = dict(zip(species,functionalForm)) 
    
    
    Temperature = modelDataObject.solution['temperature'].as_matrix()
    if physical_sens ==1:
        TempArraysForPhysicalSensitivites = modelDataObject.temps_for_physical_sens
        TempArraysForPhysicalSensitivites=TempArraysForPhysicalSensitivites[0]
        PressureArraysForPhysicalSensitivities = modelDataObject.pressure_for_physical_sens
        PressureArraysForPhysicalSensitivities=PressureArraysForPhysicalSensitivities[0]
        concentrationOfAbsorbanceObservablesForSens = modelDataObject.concentration_for_physical_sens
        concentrationOfAbsorbanceObservablesForSens = concentrationOfAbsorbanceObservablesForSens[0]

        

        dk = .01
        
        
    
    
    flat_list = [item for sublist in wavelengths for item in sublist]
    flat_list = list(set(flat_list))
    
    absorbanceSpeciesList = [[] for x in xrange(len(flat_list))]
    for i, wavelength in enumerate(flat_list):
        for j,specie in enumerate(species):
            if wavelength in SandWL[specie]:
                absorbanceSpeciesList[i].append(specie)
    
    
    d = {}
    s = {}
    ATemp= {}
    APressure = {}
    ASpecies = {}

    for i in xrange(len(absorbanceSpeciesList)):
        W = flat_list[i]
        for j in xrange(len(absorbanceSpeciesList[i])):            
            value = absorbanceSpeciesList[i][j]  
            
         
            if value == 'H2O2':
                if W == 215:
                    indx = SandWL[value].index(W)
                    ff = SandfunctionalForm[value][indx]
                    cc = SandCoupledCoefficients[value][indx]

                    if ff == 'A':
                        #print(cc[1]*Temperature)
                        epsilon = ((cc[1]*Temperature) + cc[0])*1000
                    if ff == 'B':
                        epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))*1000
                    if ff == 'C':
                        epsilon = cc[0] *1000
                    #multiplying by 1000 to convert from L to cm^3 from the epsilon given in paper 
                    #this applies if the units on epsilon are given as they are in kappl paper 
                    #must calcuate and pass in reactor volume 
                    concentration = ((np.true_divide(1,modelDataObject.solution['temperature'].as_matrix().flatten())) * (modelDataObject.solution['pressure'].as_matrix().flatten()) * (1/(8.314e6)))*modelDataObject.solution[value].as_matrix().flatten()
                    A = pathLength*(epsilon*concentration)
                    #print('this is H2o2 absorbacne')
                    #print(A)
                    if kinetic_sens == 1 :
                        
                        
                        if ff == 'C': 
                           
                            
                            absorbSens = modelDataObject.species_slice_ksens(value)*epsilon                            
                            absorbSens = absorbSens*concentration[:,np.newaxis]

                        else:

                            absorbSens = modelDataObject.species_slice_ksens(value)*epsilon[:,np.newaxis]
                            absorbSens = absorbSens*concentration[:,np.newaxis]

                        s.update({str(W)+value: absorbSens})
                    d.update({str(W) +value : A})
#editing                    
                    if physical_sens == 1:
                        for numOfPhyParams in xrange(len(physical_params)):
                            if physical_params[numOfPhyParams] == 'T':
                                indx = SandWL[value].index(W)
                                ff = SandfunctionalForm[value][indx]
                                cc = SandCoupledCoefficients[value][indx]
                                if ff == 'A':
                                    epsilon = ((cc[1]*Temperature) + cc[0])*1000
                                if ff == 'B':
                                    epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))*1000
                                if ff == 'C':
                                    epsilon = cc[0]*1000
                                concentrationT = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index('T')]))) * ((PressureArraysForPhysicalSensitivities[overallList.index('T')]) *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index('T')][value].as_matrix().flatten()
                                AT = pathLength *(epsilon*concentrationT)
                                ATemp.update({str(W)+value : AT})
                            if physical_params[numOfPhyParams] == 'P':
                                indx = SandWL[value].index(W)
                                ff = SandfunctionalForm[value][indx]
                                cc = SandCoupledCoefficients[value][indx]
                                if ff == 'A':
                                    epsilon = ((cc[1]*Temperature) + cc[0])*1000
                                if ff == 'B':
                                    epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))*1000
                                if ff == 'C':
                                    epsilon = cc[0]  *1000
                              
                                concentrationP = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index('P')]))) * ((PressureArraysForPhysicalSensitivities[overallList.index('P')]) *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index('P')][value].as_matrix().flatten()
                                AP = pathLength *(epsilon*concentrationP)
                                APressure.update({str(W) + value: AP})
                            if physical_params[numOfPhyParams] == 'X':
                                for simulationNumber, species in enumerate(physicalParamsSpecies):
                                    indx = SandWL[value].index(W)
                                    ff = SandfunctionalForm[value][indx]
                                    cc = SandCoupledCoefficients[value][indx]
                                    if ff == 'A':
                                        epsilon = ((cc[1]*Temperature) + cc[0])*1000
                                    if ff == 'B':
                                        epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))*1000
                                    if ff == 'C':
                                        epsilon = cc[0]*1000
                                    concentrationX = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index(species)]))) * (PressureArraysForPhysicalSensitivities[overallList.index(species)] *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index(species)][value].as_matrix().flatten()
                                    AX = pathLength*(epsilon * concentrationX)
                                    ASpecies.update({str(W) + value + species: AX})
                                    
                                    
                                    
                                
                    
            if value == 'HO2':
                if W == 215:
                    indx = SandWL[value].index(W)
                    ff = SandfunctionalForm[value][indx]
                    cc = SandCoupledCoefficients[value][indx]
                    if ff == 'A':
                        epsilon = ((cc[1]*Temperature) + cc[0])*1000
                    if ff == 'B':
                        epsilon = ((cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature))))))*1000
                    if ff == 'C':
                        epsilon = cc[0]  *1000

                    #this applies if the units on epsilon are given as they are in kappl paper 
                    #must calcuate and pass in reactor volume 
                    concentration = ((np.true_divide(1,modelDataObject.solution['temperature'].as_matrix().flatten())) * (modelDataObject.solution['pressure'].as_matrix().flatten()) * (1/(8.314e6)))*modelDataObject.solution[value].as_matrix().flatten()
                    A = pathLength*(epsilon*concentration)
                    #print('this is Ho2 absorbance')
                    #print(A)
                    if kinetic_sens == 1 :
                        
                        if ff == 'C':                                                  
                            absorbSens = modelDataObject.species_slice_ksens(value)*epsilon                            
                            absorbSens = absorbSens*concentration[:,np.newaxis]

                        else:
                            absorbSens = modelDataObject.species_slice_ksens(value)*epsilon[:,np.newaxis]
                            absorbSens = absorbSens*concentration[:,np.newaxis]

                        s.update({str(W)+value: absorbSens})
                    d.update({str(W) +value : A})
                    if physical_sens == 1:
                        for numOfPhyParams in xrange(len(physical_params)):
                            if physical_params[numOfPhyParams] == 'T':
                                indx = SandWL[value].index(W)
                                ff = SandfunctionalForm[value][indx]
                                cc = SandCoupledCoefficients[value][indx]
                                if ff == 'A':
                                    epsilon = ((cc[1]*Temperature) + cc[0])*1000
                                if ff == 'B':
                                    epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))*1000
                                if ff == 'C':
                                    epsilon = cc[0]*1000
                                concentrationT = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index('T')]))) * ((PressureArraysForPhysicalSensitivities[overallList.index('T')]) *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index('T')][value].as_matrix().flatten()
                                AT = pathLength *(epsilon*concentrationT)
                                ATemp.update({str(W)+value : AT})
                            if physical_params[numOfPhyParams] == 'P':
                                indx = SandWL[value].index(W)
                                ff = SandfunctionalForm[value][indx]
                                cc = SandCoupledCoefficients[value][indx]
                                if ff == 'A':
                                    epsilon = ((cc[1]*Temperature) + cc[0])*1000
                                if ff == 'B':
                                    epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))*1000
                                if ff == 'C':
                                    epsilon = cc[0]  *1000                              
                                concentrationP = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index('P')]))) * ((PressureArraysForPhysicalSensitivities[overallList.index('P')]) *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index('P')][value].as_matrix().flatten()
                                AP = pathLength *(epsilon*concentrationP)
                                APressure.update({str(W) + value: AP})
                            if physical_params[numOfPhyParams] == 'X':
                                for simulationNumber, species in enumerate(physicalParamsSpecies):
                                    indx = SandWL[value].index(W)
                                    ff = SandfunctionalForm[value][indx]
                                    cc = SandCoupledCoefficients[value][indx]
                                    if ff == 'A':
                                        epsilon = ((cc[1]*Temperature) + cc[0])*1000
                                    if ff == 'B':
                                        epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))*1000
                                    if ff == 'C':
                                        epsilon = cc[0]*1000
                                    concentrationX = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index(species)]))) * (PressureArraysForPhysicalSensitivities[overallList.index(species)] *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index(species)][value].as_matrix().flatten()
                                    AX = pathLength*(epsilon * concentrationX)
                                    ASpecies.update({str(W) + value + species: AX})                          
                    
            if value == 'HO2':
                if W == 229:
                    indx = SandWL[value].index(W)
                    ff = SandfunctionalForm[value][indx]
                    cc = SandCoupledCoefficients[value][indx]
                    if ff == 'A':
                        epsilon = (cc[1]*Temperature) + cc[0]
                    if ff == 'B':
                        epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))
                    if ff == 'C':
                        epsilon = cc[0]

                    #this applies if the units on epsilon are given as they are in kappl paper 
                    #must calcuate and pass in reactor volume 
                    concentration = ((np.true_divide(1,modelDataObject.solution['temperature'].as_matrix().flatten())) * (modelDataObject.solution['pressure'].as_matrix().flatten()) * (1/(8.314e6)))*modelDataObject.solution[value].as_matrix().flatten()                    
                    A = pathLength*(epsilon*concentration)
                    
                    if kinetic_sens == 1 :
                        
                        if ff == 'C':                                                  
                            absorbSens = modelDataObject.species_slice_ksens(value)*epsilon                            
                            absorbSens = absorbSens*concentration[:,np.newaxis]

                        else:
                            absorbSens = modelDataObject.species_slice_ksens(value)*epsilon[:,np.newaxis]
                            absorbSens = absorbSens*concentration[:,np.newaxis]

                        s.update({str(W)+value: absorbSens})
                    d.update({str(W) +value : A})                    
                    if physical_sens == 1:
                        for numOfPhyParams in xrange(len(physical_params)):
                            if physical_params[numOfPhyParams] == 'T':
                                indx = SandWL[value].index(W)
                                ff = SandfunctionalForm[value][indx]
                                cc = SandCoupledCoefficients[value][indx]
                                if ff == 'A':
                                    epsilon = (cc[1]*Temperature) + cc[0]
                                if ff == 'B':
                                    epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))
                                if ff == 'C':
                                    epsilon = cc[0]
                                concentrationT = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index('T')]))) * ((PressureArraysForPhysicalSensitivities[overallList.index('T')]) *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index('T')][value].as_matrix().flatten()
                                AT = pathLength *(epsilon*concentrationT)
                                ATemp.update({str(W)+value : AT})
                            if physical_params[numOfPhyParams] == 'P':
                                indx = SandWL[value].index(W)
                                ff = SandfunctionalForm[value][indx]
                                cc = SandCoupledCoefficients[value][indx]
                                if ff == 'A':
                                    epsilon = (cc[1]*Temperature) + cc[0]
                                if ff == 'B':
                                    epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))
                                if ff == 'C':
                                    epsilon = cc[0]                                
                                concentrationP = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index('P')]))) * ((PressureArraysForPhysicalSensitivities[overallList.index('P')]) *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index('P')][value].as_matrix().flatten()
                                AP = pathLength *(epsilon*concentrationP)
                                APressure.update({str(W) + value: AP})
                            if physical_params[numOfPhyParams] == 'X':
                                for simulationNumber, species in enumerate(physicalParamsSpecies):
                                    indx = SandWL[value].index(W)
                                    ff = SandfunctionalForm[value][indx]
                                    cc = SandCoupledCoefficients[value][indx]
                                    if ff == 'A':
                                        epsilon = (cc[1]*Temperature) + cc[0]
                                    if ff == 'B':
                                        epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))
                                    if ff == 'C':
                                        epsilon = cc[0]
                                    concentrationX = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index(species)]))) * (PressureArraysForPhysicalSensitivities[overallList.index(species)] *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index(species)][value].as_matrix().flatten()
                                    AX = pathLength*(epsilon * concentrationX)
                                    ASpecies.update({str(W) + value + species: AX})
                    
            if value == 'HO2':
                if W == 227:
                    indx = SandWL[value].index(W)
                    ff = SandfunctionalForm[value][indx]
                    cc = SandCoupledCoefficients[value][indx]
                    if ff == 'A':
                        epsilon = (cc[1]*Temperature) + cc[0]
                    if ff == 'B':
                        epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))
                    if ff == 'C':
                        epsilon = cc[0]
                    
# needed to change the divide statment to accomidate units 
#reactor volume in cm
                    concentration = ((np.true_divide(1,modelDataObject.solution['temperature'].as_matrix().flatten())) * (modelDataObject.solution['pressure'].as_matrix().flatten()) * (1/(8.314e6)))*modelDataObject.solution[value].as_matrix().flatten()
                    A = pathLength*(epsilon*concentration)
                    if kinetic_sens == 1 :
                        
                        if ff == 'C':                                                  
                            absorbSens = modelDataObject.species_slice_ksens(value)*epsilon                            
                            absorbSens = absorbSens*concentration[:,np.newaxis]

                        else:
                            absorbSens = modelDataObject.species_slice_ksens(value)*epsilon[:,np.newaxis]
                            absorbSens = absorbSens*concentration[:,np.newaxis]

                        s.update({str(W)+value: absorbSens})
                    d.update({str(W) +value : A})
                    if physical_sens == 1:
                        for numOfPhyParams in xrange(len(physical_params)):
                            if physical_params[numOfPhyParams] == 'T':
                                indx = SandWL[value].index(W)
                                ff = SandfunctionalForm[value][indx]
                                cc = SandCoupledCoefficients[value][indx]
                                if ff == 'A':
                                    epsilon = (cc[1]*Temperature) + cc[0]
                                if ff == 'B':
                                    epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))
                                if ff == 'C':
                                    epsilon = cc[0]
                                
                                
                                concentrationT = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index('T')]))) * ((PressureArraysForPhysicalSensitivities[overallList.index('T')]) *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index('T')][value].as_matrix().flatten()
                                AT = pathLength *(epsilon*concentrationT)
                                
                                ATemp.update({str(W)+value : AT})
                            if physical_params[numOfPhyParams] == 'P':
                                indx = SandWL[value].index(W)
                                ff = SandfunctionalForm[value][indx]
                                cc = SandCoupledCoefficients[value][indx]
                                if ff == 'A':
                                    epsilon = (cc[1]*Temperature) + cc[0]
                                if ff == 'B':
                                    epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))
                                if ff == 'C':
                                    epsilon = cc[0]                                
                                concentrationP = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index('P')]))) * (PressureArraysForPhysicalSensitivities[overallList.index('P')] *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index('P')][value].as_matrix().flatten()
                                AP = pathLength *(epsilon*concentrationP)
                                APressure.update({str(W) + value: AP})
                            if physical_params[numOfPhyParams] == 'X':
                                for simulationNumber, species in enumerate(physicalParamsSpecies):
                                    indx = SandWL[value].index(W)
                                    ff = SandfunctionalForm[value][indx]
                                    cc = SandCoupledCoefficients[value][indx]
                                    if ff == 'A':
                                        epsilon = (cc[1]*Temperature) + cc[0]
                                    if ff == 'B':
                                        epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))
                                    if ff == 'C':
                                        epsilon = cc[0]
                                    concentrationX = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index(species)]))) * (PressureArraysForPhysicalSensitivities[overallList.index(species)] *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index(species)][value].as_matrix().flatten()
                                    AX = pathLength*(epsilon * concentrationX)
                                    ASpecies.update({str(W) + value + species: AX})                    
            if value == 'H2O2':
                if W == 227:
                    indx = SandWL[value].index(W)
                    ff = SandfunctionalForm[value][indx]
                    cc = SandCoupledCoefficients[value][indx]
                    if ff == 'A':
                        epsilon = (cc[1]*Temperature) + cc[0]

                    if ff == 'B':
                        epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))
                    if ff == 'C':
                        epsilon = cc[0]

# needed to change the divide statment to accomidate units 
#reactor volume in cm
                    concentration = ((np.true_divide(1,modelDataObject.solution['temperature'].as_matrix().flatten())) * (modelDataObject.solution['pressure'].as_matrix().flatten()) * (1/(8.314e6)))*modelDataObject.solution[value].as_matrix().flatten()
                    A = pathLength*(epsilon*concentration)
                    if kinetic_sens == 1:
                        
                        if ff == 'C':                                                  
                            absorbSens = modelDataObject.species_slice_ksens(value)*epsilon                            
                            absorbSens = absorbSens*concentration[:,np.newaxis]

                        else:
                            absorbSens = modelDataObject.species_slice_ksens(value)*epsilon[:,np.newaxis]
                            absorbSens = absorbSens*concentration[:,np.newaxis]                    
                    
                    
                    
                        s.update({str(W)+value : absorbSens})
                    d.update({str(W) +value : A})
                    if physical_sens == 1:
                        for numOfPhyParams in xrange(len(physical_params)):
                            if physical_params[numOfPhyParams] == 'T':
                                indx = SandWL[value].index(W)
                                ff = SandfunctionalForm[value][indx]
                                cc = SandCoupledCoefficients[value][indx]
                                if ff == 'A':
                                    epsilon = (cc[1]*Temperature) + cc[0]
                                if ff == 'B':
                                    epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))
                                if ff == 'C':
                                    epsilon = cc[0]
                                
                                concentrationT = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index('T')]))) * ((PressureArraysForPhysicalSensitivities[overallList.index('T')]) *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index('T')][value].as_matrix().flatten()
                                AT = pathLength *(epsilon*concentrationT)
                                ATemp.update({str(W)+value : AT})
                                
                            if physical_params[numOfPhyParams] == 'P':
                                indx = SandWL[value].index(W)
                                ff = SandfunctionalForm[value][indx]
                                cc = SandCoupledCoefficients[value][indx]
                                if ff == 'A':
                                    epsilon = (cc[1]*Temperature) + cc[0]
                                if ff == 'B':
                                    epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))
                                if ff == 'C':
                                    epsilon = cc[0]                                
                                concentrationP = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index('P')]))) * ((PressureArraysForPhysicalSensitivities[overallList.index('P')]) *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index('P')][value].as_matrix().flatten()
                                AP = pathLength *(epsilon*concentrationP)
                                APressure.update({str(W) + value: AP})
                            if physical_params[numOfPhyParams] == 'X':
                                for simulationNumber, species in enumerate(physicalParamsSpecies):
                                    indx = SandWL[value].index(W)
                                    ff = SandfunctionalForm[value][indx]
                                    cc = SandCoupledCoefficients[value][indx]
                                    if ff == 'A':
                                        epsilon = (cc[1]*Temperature) + cc[0]
                                    if ff == 'B':
                                        epsilon = (cc[0]*(1-(np.exp(np.true_divide(cc[1],Temperature)))))
                                    if ff == 'C':
                                        epsilon = cc[0]
                                    concentrationX = ((np.true_divide(1,TempArraysForPhysicalSensitivites[overallList.index(species)]))) * (PressureArraysForPhysicalSensitivities[overallList.index(species)] *(1/(8.314e6))) * concentrationOfAbsorbanceObservablesForSens[overallList.index(species)][value].as_matrix().flatten()
                                    AX = pathLength*(epsilon * concentrationX)
                                    ASpecies.update({str(W) + value + species: AX})
                        
            
           
            
    aslWithW = [[] for x in xrange(len(flat_list))]
    for i, wave in enumerate(flat_list):
        for j, specie in enumerate(absorbanceSpeciesList[i]):
            aslWithW[i].append(str(wave)+specie)
    #not right         
    aslWithWWithSpecies = []
    for x, specie in enumerate(physicalParamsSpecies):
        for i in xrange(len(aslWithW)):
            tempx=[]
            for y, combo in enumerate(aslWithW[i]):
                tempx.append(combo+specie)
            aslWithWWithSpecies.append(tempx)

    absorbances = [[] for x in xrange(len(flat_list))]
    for x in xrange(len(aslWithW)):
        for i, indx in enumerate(aslWithW[x]):
            absorbances[x].append(d[indx])
            
    if physical_sens ==1: 
        
        absorbancesT = [[] for x in xrange(len(flat_list))]
        for x in xrange(len(aslWithW)):
            for i, indx in enumerate(aslWithW[x]):
                absorbancesT[x].append(ATemp[indx])

        
                
        absorbancesP = [[] for x in xrange(len(flat_list))]
        for x in xrange(len(aslWithW)):
            for i, indx in enumerate(aslWithW[x]):
                absorbancesP[x].append(APressure[indx])  
               
        absorbancesX = [[] for x in xrange(len(flat_list)*len(physicalParamsSpecies))]        
        for x in xrange(len(aslWithWWithSpecies)):
            for y, indx in enumerate(aslWithWWithSpecies[x]):
                    absorbancesX[x].append(ASpecies[indx])

      
    sensitivities = [[] for x in xrange(len(flat_list))]
    for x in xrange(len(aslWithW)):
        for i,indx in enumerate(aslWithW[x]):
            sensitivities[x].append(s[indx])

    
    
           
    absorbanceSum=[]
    for i in xrange(len(flat_list)):
        absorbanceSum.append(sum(absorbances[i]))
    if physical_sens ==1:    
        absorbanceTSum=[]
        for i in xrange(len(flat_list)):
            absorbanceTSum.append(sum(absorbancesT[i]))
            
        absorbancePSum =[]
        for i in xrange(len(flat_list)):
            absorbancePSum.append(sum(absorbancesP[i]))
            
        absorbanceXSum = [[] for x in xrange(len(flat_list)*len(physicalParamsSpecies))]
        for i in xrange(len(flat_list)*len(physicalParamsSpecies)):
            absorbanceXSum[i].append(sum(absorbancesX[i]))

    
        
    absorbanceSensSum = []
    for i in xrange(len(flat_list)):
        absorbanceSensSum.append(sum(sensitivities[i]))

        
    LoverA = [np.true_divide(pathLength, absorbanceSum[WLs]) for WLs in xrange(len(absorbanceSum)) ]

    #x=0
    absorbanceSensSum = [absorbanceSensSum[x]*LoverA[x][:,np.newaxis] for x in xrange(len(absorbanceSum))]

    ##rearange absorbance X sum
    if physical_sens ==1:
        absorbanceXSumFlattened = [item for sublist in absorbanceXSum for item in sublist]
        def slice_per(source,step):
            return [source[i::step] for i in range(step)]
        absorbanceXSumReArranged = slice_per(absorbanceXSumFlattened,len(flat_list)) 
    
    print(physicalParamsSpecies)
    print('this is pp')
    indexForX = []        
    for WL, physicalParam in itertools.product(flat_list,physicalParamsSpecies):
        indexForX.append(str(WL)+physicalParam)
    print(indexForX)
    #indexForX = [indexForX[i:i+3] for i in range(0, len(indexForX), len(physicalParam))]
    indexForX = [indexForX]             
    print('this is the index for x')
    print(indexForX)


    
            
    WLandAbsSum = dict(zip(flat_list,absorbanceSum))
    WLandSensSum = dict(zip(flat_list,absorbanceSensSum))
    absorbanceSensAsum = absorbanceSensSum
    absorbanceSensNsum = [absorbanceSensSum[x]* (np.log(modelDataObject.solution['temperature'].as_matrix().flatten()))[:,np.newaxis] for x in xrange(len(absorbanceSensSum))]
    absorbanceSensEasum = [absorbanceSensSum[x]* 1/ct.gas_constant*(np.true_divide(-1,modelDataObject.solution['temperature'].as_matrix().flatten()))[:,np.newaxis] for x in xrange(len(absorbanceSensSum))]
        
    WLandSensAsum = dict(zip(flat_list,absorbanceSensAsum))
    WLandSensNsum = dict(zip(flat_list,absorbanceSensNsum))
    WLandSensEasum = dict(zip(flat_list,absorbanceSensEasum))
    
    
    
   
    
    if physical_sens ==1:
        WLandAbsTsum = dict(zip(flat_list,absorbanceTSum))
        WLandAbsPsum = dict(zip(flat_list,absorbancePSum))
        WLandAbsXsum = dict(zip(flat_list,absorbanceXSumReArranged))
        
        
    
    
    #mapping absorbance to A, n and Ea



    
    #start interpolation editing  
            
            
    absorbanceExperimentalData = [pd.read_csv(csv) for csv in absorbanceCsvFiles]
    absorbanceExperimentalData = [absorbanceExperimentalData[x].dropna(how = 'any') for x in xrange(len(absorbanceExperimentalData))]
    absorbanceExperimentalData = [absorbanceExperimentalData[x].apply(pd.to_numeric, errors = 'coerce').dropna() for x in xrange(len(absorbanceExperimentalData))]
    absorbanceExperimentalDataList = [absorbanceExperimentalData[x].as_matrix() for x in xrange(len(absorbanceExperimentalData))]                                    
    #time = modelDataObject.solution['time']
    time = [absorbanceExperimentalData[x].ix[:,0].values for x in xrange(len(absorbanceExperimentalData))]

    



    #interpolatedAbsorbanceData = [np.interp(time,absorbanceExperimentalData[x].ix[:,0].values,absorbanceExperimentalData[x].ix[:,1].values) for x in xrange(len(absorbanceExperimentalData))]
    interpolatedAbsorbanceData = [np.interp(time[x],modelDataObject.solution['time'].values,WLandAbsSum[waveLength]) for x,waveLength in enumerate(flat_list)]

    interpolatedAbsorbanceData = [pd.DataFrame(interpolatedAbsorbanceData[x]) for x in xrange(len(interpolatedAbsorbanceData))]
    
    interpolatedAbsorbanceData = [interpolatedAbsorbanceData[x].as_matrix().flatten() for x in xrange(len(interpolatedAbsorbanceData))] 
    interpolatedAbsorbanceDataDict = dict(zip(flat_list,interpolatedAbsorbanceData))
    modelDataObject.add_interpolated_simulation_data_absorbance(interpolatedAbsorbanceDataDict)
           
    
    diffOfPredictAndExpA = []
            
    for x,WL in enumerate(absorbanceCsvWavelengths):
        diffOfPredictAndExpA.append(np.subtract(np.log(absorbanceExperimentalData[x].ix[:,1].values),np.log(interpolatedAbsorbanceData[x])))
        
    diffOfPredictAndExpA = [pd.DataFrame(diffOfPredictAndExpA[x]) for x in xrange(len(diffOfPredictAndExpA))]
    diffOfPredictAndExpA = [diffOfPredictAndExpA[x].as_matrix().flatten() for x in xrange(len(diffOfPredictAndExpA))]
    diffOfPredictAndExpA = [diffOfPredictAndExpA[x].reshape(int(np.shape(diffOfPredictAndExpA[x])[0]),1)]

    
    if len(diffOfPredictAndExpA) > 1:                 
        diffOfPredictAndExpA = np.column_stack((diffOfPredictAndExpA))
    else: 
        diffOfPredictAndExpA = np.array(diffOfPredictAndExpA[0])
        
        

                   
    interpolatedAbsorbanceData = [pd.DataFrame(interpolatedAbsorbanceData[x]) for x in xrange(len(interpolatedAbsorbanceData))]
    interpolatedAbsorbanceData = [interpolatedAbsorbanceData[x].as_matrix().flatten() for x in xrange(len(interpolatedAbsorbanceData))]
    interpolatedAbsorbanceData = [interpolatedAbsorbanceData[x].reshape(int(np.shape(interpolatedAbsorbanceData[x])[0]),1)]
    
    
    if len(interpolatedAbsorbanceData) > 1:
        interpolatedAbsorbanceData = np.column_stack((interpolatedAbsorbanceData))
    else:
        interpolatedAbsorbanceData = np.array(interpolatedAbsorbanceData[0])
        
   #editing here ...  for doing absorbance sensitiviteis 

     
    interpolatedSensitivitysAA = []
    for x,sensitivity in enumerate(flat_list):
        temp = []
        for column in xrange(int(np.shape(WLandSensAsum[sensitivity])[1])):
            interpolatedColumn = np.interp(time[x],modelDataObject.solution['time'].values,WLandSensAsum[sensitivity][:,column])
            temp.append(interpolatedColumn.reshape(int(np.shape(interpolatedColumn)[0]),1))
        temp = np.column_stack(temp)
        interpolatedSensitivitysAA.append(temp)

    
    interpolatedSensitivitysAn = []
    for x,sensitivity in enumerate(flat_list):
        temp2 = []
        for column in xrange(int(np.shape(WLandSensNsum[sensitivity])[1])):
            interpolatedColumn = np.interp(time[x],modelDataObject.solution['time'].values,WLandSensNsum[sensitivity][:,column])
            temp2.append(interpolatedColumn.reshape(int(np.shape(interpolatedColumn)[0]),1))
        temp2 = np.column_stack(temp2)
        interpolatedSensitivitysAn.append(temp2)               
                
        
    interpolatedSensitivitysAEa = []
    for x,sensitivity in enumerate(flat_list):
        temp3 = []
        for column in xrange(int(np.shape(WLandSensEasum[sensitivity])[1])):
            interpolatedColumn = np.interp(time[x],modelDataObject.solution['time'].values,WLandSensEasum[sensitivity][:,column])
            temp3.append(interpolatedColumn.reshape(int(np.shape(interpolatedColumn)[0]),1))
        temp3 = np.column_stack(temp3)
        interpolatedSensitivitysAEa.append(temp3) 
        
        
        
    InterpolatedWLandSensAsum  = dict(zip(flat_list,interpolatedSensitivitysAA))
    WLandSensAsum = dict(zip(flat_list,interpolatedSensitivitysAA))
    WLandSensNsum = dict(zip(flat_list,interpolatedSensitivitysAn))
    WLandSensEasum = dict(zip(flat_list,interpolatedSensitivitysAEa))
    
        
        
            
        #zip this back up so I don't have to change that much about the S matrix??
        
        #now going to need to write this to the simulations model data object
        

        
     # interpolated sensitivity    
        
    #calcualte absorbance sensitivity for T , P , X
    if physical_sens ==1:
        AbsorbanceSensForT = [[] for x in xrange(len(absorbanceCsvWavelengths))]
        for x,WL in enumerate(absorbanceCsvWavelengths):
            #AbsorbanceSensForT[x].append(np.true_divide(np.subtract(np.log10(WLandAbsSum[WL]),np.log10(WLandAbsTsum[WL])),np.log10(dk)))
            AbsorbanceSensForT[x].append(np.true_divide(np.subtract(np.log(WLandAbsTsum[WL]),np.log(WLandAbsSum[WL])),dk))
            
            
        AbsorbanceSensForP = [[] for x in xrange(len(absorbanceCsvWavelengths))]    
        for x,WL in enumerate(absorbanceCsvWavelengths):
            #AbsorbanceSensForP[x].append(np.true_divide(np.subtract(np.log10(WLandAbsSum[WL]),np.log10(WLandAbsPsum[WL])),np.log10(dk)))
            AbsorbanceSensForP[x].append(np.true_divide(np.subtract(np.log(WLandAbsPsum[WL]),np.log(WLandAbsSum[WL])),dk))
    
        AbsorbanceSensForX = [[] for x in xrange(len(absorbanceCsvWavelengths))]
        for x, WL in enumerate(flat_list):
            
            for y, combo in enumerate(indexForX[x]):
                print(y,combo)
                AbsorbanceSensForX[x].append(np.true_divide(np.subtract(np.log(WLandAbsXsum[WL][y]),np.log(WLandAbsSum[WL])),dk))
                
        AbsorbanceSensForXCombined = [[] for x in xrange(len(absorbanceCsvWavelengths))]
        for x in xrange(len(absorbanceCsvWavelengths)):
            AbsorbanceSensForXCombined[x].append(np.column_stack((AbsorbanceSensForX[x])))
        
            
        print(np.shape(AbsorbanceSensForXCombined[0][0]))
    

        AbsorbancePSensCombinedByWL = [[] for x in xrange(len(absorbanceCsvWavelengths))]
        for WL in xrange(len(absorbanceCsvWavelengths)):    
            AbsorbancePSensCombinedByWL[WL].append(np.column_stack((np.transpose(AbsorbanceSensForT[x]), np.transpose(AbsorbanceSensForP[x]), AbsorbanceSensForXCombined[x][0])))

        print('abs combined')
        print(np.shape(AbsorbancePSensCombinedByWL[0][0]))    
        interpolatedPhysicalSensitivities = [[] for x in xrange(len(absorbanceCsvWavelengths))]
        for x,sensitivity in enumerate(flat_list):
            temp4 = []
            for column in xrange(int(np.shape(AbsorbancePSensCombinedByWL[x][0])[1])):
                interpolatedColumn = np.interp(time[x],modelDataObject.solution['time'].values,AbsorbancePSensCombinedByWL[x][0][:,column])
                temp4.append(interpolatedColumn.reshape(int(np.shape(interpolatedColumn)[0]),1))                
            temp4 = np.column_stack(temp4)
            interpolatedPhysicalSensitivities[x].append(temp4)
        
            


    
            
    print('this ist he array shape being passed')        
    print(np.shape(interpolatedPhysicalSensitivities[0][0]))       
    AExpIndex = [modelDataObject.solution['time'].as_matrix(),absorbanceCsvWavelengths]
    #do we want to return individual species that make up the absorbance curve?
    #changd this to return absorbanceexperimentaldata list 
    modelDataObject.add_exp_shocktube_absorbance(absorbanceExperimentalDataList,diffOfPredictAndExpA,absorbanceProfiles = WLandAbsSum, absorbanceSens=WLandSensSum ,index=AExpIndex)
    
    modelDataObject.absorbance_sens_mappings( WLandSensAsum , WLandSensNsum, WLandSensEasum) 
    if physical_sens ==1:
        modelDataObject.add_exp_shocktube_absorbance(absorbanceExperimentalDataList,diffOfPredictAndExpA,absorbanceProfiles = WLandAbsSum, absorbanceSens=WLandSensSum ,index=AExpIndex, absorbancePSens=interpolatedPhysicalSensitivities)
        
       # modelDataObject.absorbancePSens=AbsorbancePSensCombinedByWL
        
        
    speciesNames.append(addBackin)
    return modelDataObject
