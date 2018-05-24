# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 11:10:45 2017

@author: Mark Barbet
"""

#Code is for a free flame simulation with and without sensitivities

from __future__ import print_function
from __future__ import division

import numpy as np
import cantera as ct
import pandas as pd
import simulations as sim


def free_flame2(gas,width,reactants):
    #Inlet Temperature in Kelvin and Inlet Pressure in Pascals
    #In this case we are setting the inlet T and P to room temperature conditions
    To = gas.T
    Po = gas.P

    #Define the gas-mixutre and kinetics
    #In this case, we are choosing a GRI3.0 gas
    gas = ct.Solution('gri30.cti')

    # Create a stoichiometric CH4/Air premixed mixture 
    gas.set_equivalence_ratio(1.0, 'CH4', {'O2':1.0, 'N2':3.76})
    gas.TP = To, Po
    
    # Domain width in metres
    width = 0.114

    # Create the flame object
    flame = ct.FreeFlame(gas, width=width)

    # Define tolerances for the solver
    flame.set_refine_criteria(ratio=2, slope=0.1, curve=0.1)

    # Define logging level
    loglevel = 1
    flame.solve(loglevel=loglevel, auto=True)
    Su0 = flame.u[0]
    print("Flame Speed is: {:.2f} cm/s".format(Su0*100))

    #Note that the variable Su0 will also be used downsteam in the sensitivity analysis
    
    
def free_flame(phi,fuel,reactants,gas,width,data=pd.DataFrame(columns=['z','T']),kinetic_sens=0,physical_sens=0,observables=[],physical_params=['T','P'],energycon=False,soret=True,flamespeed_sens=0):
    #when energycon is off treat flame as burner stabilized, with known T-profile
    gas.set_equivalence_ratio(phi, fuel, reactants)
    simtype = 'free flame'
    baseConditions=gas.TPX
    tol_ss = [1.0e-5, 1.0e-13]  # [rtol atol] for steady-state problem
    tol_ts = [1.0e-4, 1.0e-10]  # [rtol atol] for time stepping
    loglevel = 1  # amount of diagnostic output (0 to 5)
   
    #f = ct.FreeFlame(gas, grid=np.linspace(0,width,num=20))
    f = ct.FreeFlame(gas,width=width)
    #f.burner.mdot = mdot
    
    #f.inlet.T = baseConditions[0]
    #f.inlet.X = reactants
    #f.set_initial_guess()
    
    #f.set_initial_guess()
    # read temperature vs. position data from a file.
    # The file is assumed to have one z, T pair per line, separated by a comma.
    if data.empty==False and energycon==False:
        zloc=data['z']
        tvalues=data['T']
        #zloc, tvalues = np.genfromtxt('C:\\Users\\HP USER\\Desktop\\kevin\\harrington_78torr_Tprofile_chemkin.txt', delimiter=',', comments='#').T
        #zloc, tvalues = np.genfromtxt('C:\\Users\\HP USER\\Desktop\\kevin\\t_profile_38torr.txt', delimiter=',', comments='#').T    
        zloc /= max(zloc)

        #print(tvalues)

        f.flame.set_fixed_temp_profile(zloc, tvalues)  #sets a fixed temperature profile for the flame simulation.  May come from a measurement. Requires no energy conservation
    elif data.empty==False and energycon!='off':
        raise Exception('User has supplied fixed temperature dataset but energy conservation is not off.  Remove dataset or turn energy conservation off')
        
        
    f.flame.set_steady_tolerances(default=tol_ss)  #Set steady tolerances
    f.flame.set_transient_tolerances(default=tol_ts) #Set transient tolerances
    #f.show_solution()

    f.energy_enabled = energycon  #This must be set to false for a burner stabilized flame with known T-profile
    #f.max_time_step_count=1500
    
    
    f.transport_model = 'Multi'   #Sets to multicomponent transport for simulation.  Needs to be set this way to use Soret effect
    f.set_max_jac_age(10, 10)       #Age limits on Jacobian-leave as is for best results
    f.solve(loglevel, refine_grid=False)  #Solve for initial estimate without grid refinement
    #f.transport_model = 'Mix'
    f.soret_enabled = soret          #Enable Soret effect.  Remember transport must be set to multi.  Mix causes failure
    
    f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.1)  #Establishes refinement criteria for grid
    #f.solve(loglevel,refine_grid=True)
    #print('mixture-averaged flamespeed = ', f.u[0])
    f.transport_model = 'Multi'
             #This block solves problem again with grid refinement on
    f.solve(loglevel, refine_grid=True)
    forward_rates=f.forward_rate_constants
    net_rates=f.net_rates_of_progress
    reverse_rates=f.reverse_rate_constants
    #f.solve(loglevel,auto=True)
    #f.show_solution()
    #print('multicomponent flamespeed = ', f.u[0])
    Su0 = f.u[0]

    
    #solution = f
    ##Begin section to calculate sensitivities
    dk = 0.010
    solution = f.X
    if kinetic_sens==1 and bool(observables):
        
        #Calculate kinetic sensitivities
        sensIndex = [f.grid.tolist(),gas.reaction_equations(),observables]
        
        S = np.zeros((len(f.grid),gas.n_reactions,len(observables)))
        #print(solution.X[solution.flame.component_index(observables[0])-4,len(f.grid)-1])
        #a=solution.X[solution.flame.component_index(observables[0])-4,len(f.grid)-1]
        for m in range(gas.n_reactions):
            gas.set_multiplier(1.0)
            gas.set_multiplier(1+dk,m)
            f.solve(loglevel=1,refine_grid=False)
            for i in np.arange(len(observables)):
                for k in np.arange(len(f.grid)):                    
                    S[k,m,i]=np.log10(solution[f.flame.component_index(observables[i])-4,k])-np.log10(f.X[f.flame.component_index(observables[i])-4,k])
                    #print(solution.X[solution.flame.component_index(observables[i])-4,k])
                    #print(f.X[f.flame.component_index(observables[i])-4,k])
                    S[k,m,i]=np.divide(S[k,m,i],np.log10(dk))             
                    
    
#    if flamespeed_sens==1:
#        fsens = pd.DataFrame(data=[], index=gas.reaction_equations(range(gas.n_reactions)))
#        dk = 1e-2   
#        fsens["baseCase"] = ""  
#        for m in range(gas.n_reactions):
#            gas.set_multiplier(1.0) # reset all multipliers                                                                     
#            gas.set_multiplier(1+dk, m) # perturb reaction m   
#    
#            # Always force loglevel=0 for this
#            # Make sure the grid is not refined, otherwise it won't strictly 
#            # be a small perturbation analysis
#            f.solve(loglevel=0, refine_grid=False)
#    
#            # The new flame speed
#            Su = f.u[0]
#    
#            fsens["baseCase"][m] = (Su-Su0)/(Su0*dk)
    if flamespeed_sens==1:    
        fsens = f.get_flame_speed_reaction_sensitivities()

        # This step is essential, otherwise the mechanism will have been altered
        gas.set_multiplier(1.0)
    
    if physical_sens==1 and bool(observables):
        #Calculate physical sensitivities
        gas.set_multiplier(1.0)
        
        psensIndex = [f.grid.tolist(),physical_params,observables]
        pS = np.zeros((len(f.grid),len(physical_params),len(observables)))
        for m in range(len(physical_params)):
            gas.TPX=baseConditions
            if physical_params[m]=='T':
                gas.TPX=baseConditions[0]+dk,baseConditions[1],baseConditions[2]
            elif physical_params[m]=='P':
                gas.TPX=baseConditions[0],baseConditions[1]+dk,baseConditions[2]
            f.solve(loglevel=1,refine_grid=False)
            for i in np.arange(len(observables)):
                for k in np.arange(len(f.grid)):
                    pS[k,m,i] =np.log10(solution[f.flame.component_index(observables[i])-4,k])-np.log10(f.X[f.flame.component_index(observables[i])-4,k])
                    pS[k,m,i] = np.divide(pS[k,m,i],np.log10(dk))
                    
                    
                    
    
    elif kinetic_sens==1 and bool(observables)==False:
        raise Exception('Please supply a list of observables in order to run kinetic sensitivity analysis')
    elif physical_sens==1 and bool(observables)==False:
        raise Exception('Please supply a list of observables in order to run physical sensitivity analysis')
    gas.set_multiplier(1.0)
    gas.TP = baseConditions[0],baseConditions[1]  
    f.solve(loglevel=1,refine_grid=False)
    solution=pd.DataFrame(columns=f.flame.component_names)
    for i in f.flame.component_names:
        solution[i]=f.solution(i)
    
       
    if kinetic_sens==1 and bool(observables) and physical_sens!=1:
        results = sim.model_data(simtype,kinetic_sens=S,Solution=solution,Index=sensIndex)        
        
    elif physical_sens==1 and bool(observables) and kinetic_sens!=1:
        results = sim.model_data(simtype,Solution=solution,pIndex=psensIndex,physical_sens=pS)
        
    elif kinetic_sens==1 and physical_sens==1 and bool(observables):
        results = sim.model_data(simtype,Solution=solution,pIndex=psensIndex,Index=sensIndex,physical_sens=pS,kinetic_sens=S)
        
    elif kinetic_sens!=1 and physical_sens!=1:
        Index = [f.grid.tolist(),gas.reaction_equations()] 
        results = sim.model_data(simtype,Solution=solution,Index=Index)
        
    else:
        print('Something went wrong with the parameters given.  Kinetic and physical sens may be set to either 0 or 1, and require a list of observables')
    if flamespeed_sens==1:
        results.assign_flamespeed_sens(fsens,gas.reaction_equations())
    results.assign_phi(phi)
    results.add_fuel(fuel)
    results.add_forward_rates(forward_rates)
    results.add_net_rates_of_progress(net_rates)
    results.add_reverse_rates(reverse_rates)
    return results
            