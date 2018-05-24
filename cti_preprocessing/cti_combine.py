# -*- coding: utf-8 -*-
"""
Created on Thu Aug 03 14:32:24 2017

@author: Mark Barbet
"""
"""Active parameter CTI writer.  Function takes a subset of reactions from an already modified cti file and
writes them to internal memory.  It then reads an input from a portion of the code dealing with master equation simulation 
and adds those reactions to create a complete internal mechanism
"""

import numpy as np
import cantera as ct
import soln2cti as ctiw


def cti_write(x={},original_cti='',master_rxns='',master_index=[]):
    if not original_cti:
        raise Exception('Please provide a name for the original mechanism file and try again.')
    if not master_rxns and np.any(master_index):
        raise Exception('Please provide a mechanism file for reactions analysed with master equation or leave master_index empty')
    if master_rxns and not np.any(master_index):
        raise Exception('Please provide master_index, a non-empty list of reaction numbers from original file which are analysed with master equation.')
        
    if not master_rxns and not master_index:
        master_index=np.ones(ct.Solution(original_cti).n_reactions,dtype=bool)
    elif master_rxns and np.any(master_index):
        temp=np.ones(ct.Solution(original_cti).n_reactions,dtype=bool)
        for j in np.arange(len(master_index)):
            
            temp[master_index[j]-1]=False        
        master_index=temp
    lineList=[]
    with open(original_cti) as f:
        lineList=f.readlines()        
    done=False
    count=0
    while not done or count<len(lineList):
        if 'Reaction data' in lineList[count] or 'Reaction Data' in lineList[count] or 'reaction data' in lineList[count]:
            done=True
            lineList=lineList[0:count-1]
        else:count+=1
    with open('tempcti.cti','w') as p:
        p.writelines(lineList)
        
    NewModel=ct.Solution('tempcti.cti')
    original_mechanism=ct.Solution(original_cti)
    master_count=0
    if master_rxns:
        master_reactions=ct.Solution(master_rxns)
    for i in np.arange(original_mechanism.n_reactions):
        if master_index[i]==True:
            NewModel.add_reaction(original_mechanism.reaction(i))
        elif master_index[i]==False:
            NewModel.add_reaction(master_reactions.reaction(master_count))
            master_count+=1   
    
    if x=={}:
        for j in np.arange(original_mechanism.n_reactions):
           
           if master_index[j]:
               if 'ThreeBodyReaction' in str(type(original_mechanism.reaction(j))):
                   NewModel.reaction(j).rate=original_mechanism.reaction(j).rate
               elif 'ElementaryReaction' in str(type(original_mechanism.reaction(j))):
                   NewModel.reaction(j).rate=original_mechanism.reaction(j).rate
               elif 'FalloffReaction' in str(type(original_mechanism.reaction(j))):
                   NewModel.reaction(j).high_rate=original_mechanism.reaction(j).high_rate
                   NewModel.reaction(j).low_rate=original_mechanism.reaction(j).low_rate
                   
                   if original_mechanism.reaction(j).falloff.type=='Troe':
                       NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                   if original_mechanism.reaction(j).falloff.type=='Sri':
                       NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
               elif 'ChemicallyActivatedReaction' in str(type(original_mechanism.reaction(j))):
                   NewModel.reaction(j).high_rate=original_mechanism.reaction(j).high_rate
                   NewModel.reaction(j).low_rate=original_mechanism.reaction(j).low_rate
                   if original_mechanism.reaction(j).falloff.type=='Troe':
                       NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                   if original_mechanism.reaction(j).falloff.type=='Sri':
                       NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
               elif 'PlogReaction' in str(type(original_mechanism.reaction(j))):
                   NewModel.reaction(j).rates=original_mechanism.reaction(j).rates
               elif 'ChebyshevReaction' in str(type(original_mechanism.reaction(j))):
                   NewModel.reaction(j).set_parameters(original_mechanism.reaction(j).Tmin,original_mechanism.reaction(j).Tmax,original_mechanism.reaction(j).Pmin,original_mechanism.reaction(j).Pmax,original_mechanism.reaction(j).coeffs)
    if x!={}:
        for j in np.arange(original_mechanism.n_reactions):
           if master_index[j]:
               try:
                   if 'ThreeBodyReaction' in str(type(original_mechanism.reaction(j))):
                       A=original_mechanism.reaction(j).rate.pre_exponential_factor
                       n=original_mechanism.reaction(j).rate.temperature_exponent
                       Ea=original_mechanism.reaction(j).rate.activation_energy
                       NewModel.reaction(j).rate=ct.Arrhenius(A*np.exp(x['r'+str(j)]['A']),n+x['r'+str(j)]['n'],Ea+x['r'+str(j)]['Ea'])
                   elif 'ElementaryReaction' in str(type(original_mechanism.reaction(j))):
                       A=original_mechanism.reaction(j).rate.pre_exponential_factor
                       n=original_mechanism.reaction(j).rate.temperature_exponent
                       Ea=original_mechanism.reaction(j).rate.activation_energy
                       NewModel.reaction(j).rate=ct.Arrhenius(A*np.exp(x['r'+str(j)]['A']),n+x['r'+str(j)]['n'],Ea+x['r'+str(j)]['Ea'])
                   elif 'FalloffReaction' in str(type(original_mechanism.reaction(j))):
                       A=original_mechanism.reaction(j).high_rate.pre_exponential_factor
                       n=original_mechanism.reaction(j).high_rate.temperature_exponent
                       Ea=original_mechanism.reaction(j).high_rate.activation_energy
                       NewModel.reaction(j).high_rate=ct.Arrhenius(A*np.exp(x['r'+str(j)]['A']),n+x['r'+str(j)]['n'],Ea+x['r'+str(j)]['Ea'])
                       A=original_mechanism.reaction(j).low_rate.pre_exponential_factor
                       n=original_mechanism.reaction(j).low_rate.temperature_exponent
                       Ea=original_mechanism.reaction(j).low_rate.activation_energy
                       NewModel.reaction(j).low_rate=ct.Arrhenius(A*np.exp(x['r'+str(j)]['A']),n+x['r'+str(j)]['n'],Ea+x['r'+str(j)]['Ea'])
                       if original_mechanism.reaction(j).falloff.type=='Troe':
                           NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                       if original_mechanism.reaction(j).falloff.type=='Sri':
                           NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                   elif 'ChemicallyActivatedReaction' in str(type(original_mechanism.reaction(j))):
                       A=original_mechanism.reaction(j).high_rate.pre_exponential_factor
                       n=original_mechanism.reaction(j).high_rate.temperature_exponent
                       Ea=original_mechanism.reaction(j).high_rate.activation_energy
                       NewModel.reaction(j).high_rate=ct.Arrhenius(A*np.exp(x['r'+str(j)]['A']),n+x['r'+str(j)]['n'],Ea+x['r'+str(j)]['Ea'])
                       A=original_mechanism.reaction(j).low_rate.pre_exponential_factor
                       n=original_mechanism.reaction(j).low_rate.temperature_exponent
                       Ea=original_mechanism.reaction(j).low_rate.activation_energy
                       NewModel.reaction(j).low_rate=ct.Arrhenius(A*np.exp(x['r'+str(j)]['A']),n+x['r'+str(j)]['n'],Ea+x['r'+str(j)]['Ea'])
                       if original_mechanism.reaction(j).falloff.type=='Troe':
                           NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                       if original_mechanism.reaction(j).falloff.type=='Sri':
                           NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                   elif 'PlogReaction' in str(type(original_mechanism.reaction(j))):
                       NewModel.reaction(j).rates=original_mechanism.reaction(j).rates
                   elif 'ChebyshevReaction' in str(type(original_mechanism.reaction(j))):
                       NewModel.reaction(j).set_parameters(original_mechanism.reaction(j).Tmin,original_mechanism.reaction(j).Tmax,original_mechanism.reaction(j).Pmin,original_mechanism.reaction(j).Pmax,original_mechanism.reaction(j).coeffs)
               except:
                   if 'ThreeBodyReaction' in str(type(original_mechanism.reaction(j))):
                       NewModel.reaction(j).rate=original_mechanism.reaction(j).rate
                   elif 'ElementaryReaction' in str(type(original_mechanism.reaction(j))):
                       NewModel.reaction(j).rate=original_mechanism.reaction(j).rate
                   elif 'FalloffReaction' in str(type(original_mechanism.reaction(j))):
                       NewModel.reaction(j).high_rate=original_mechanism.reaction(j).high_rate
                       NewModel.reaction(j).low_rate=original_mechanism.reaction(j).low_rate
                       if original_mechanism.reaction(j).falloff.type=='Troe':
                           NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                       if original_mechanism.reaction(j).falloff.type=='Sri':
                           NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                   elif 'ChemicallyActivatedReaction' in str(type(original_mechanism.reaction(j))):
                      NewModel.reaction(j).high_rate=original_mechanism.reaction(j).high_rate
                      NewModel.reaction(j).low_rate=original_mechanism.reaction(j).low_rate
                      if original_mechanism.reaction(j).falloff.type=='Troe':
                          NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                      if original_mechanism.reaction(j).falloff.type=='Sri':
                          NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                   elif 'PlogReaction' in str(type(original_mechanism.reaction(j))):
                      NewModel.reaction(j).rates=original_mechanism.reaction(j).rates
                   elif 'ChebyshevReaction' in str(type(original_mechanism.reaction(j))):
                      NewModel.reaction(j).set_parameters(original_mechanism.reaction(j).Tmin,original_mechanism.reaction(j).Tmax,original_mechanism.reaction(j).Pmin,original_mechanism.reaction(j).Pmax,original_mechanism.reaction(j).coeffs)
                   
               
    
    new_file=ctiw.write(NewModel)
    return new_file
def cti_write2(x={},original_cti='',master_rxns='',master_index=[]):
    if not original_cti:
        raise Exception('Please provide a name for the original mechanism file and try again.')
    if not master_rxns and np.any(master_index):
        raise Exception('Please provide a mechanism file for reactions analysed with master equation or leave master_index empty')
    if master_rxns and not np.any(master_index):
        raise Exception('Please provide master_index, a non-empty list of reaction numbers from original file which are analysed with master equation.')
        
    if not master_rxns and not master_index:
        master_index=np.ones(ct.Solution(original_cti).n_reactions,dtype=bool)
    elif master_rxns and np.any(master_index):
        temp=np.ones(ct.Solution(original_cti).n_reactions,dtype=bool)
        for j in np.arange(len(master_index)):
            
            temp[master_index[j]-1]=False        
        master_index=temp
    lineList=[]
    with open(original_cti) as f:
        lineList=f.readlines()        
    done=False
    count=0
    while not done or count<len(lineList):
        if 'Reaction data' in lineList[count] or 'Reaction Data' in lineList[count] or 'reaction data' in lineList[count]:
            done=True
            lineList=lineList[0:count-1]
        else:count+=1
    with open('tempcti.cti','w') as p:
        p.writelines(lineList)
    
    NewModel=ct.Solution('tempcti.cti')
    original_mechanism=ct.Solution(original_cti)
    original_rxn_count=0
    master_rxn_eqs=[]
    if master_rxns:
        with open(master_rxns) as f:
            reactionsList=f.readlines()
        lineList=lineList+reactionsList
        with open('masterTemp.cti') as f:
            f.writelines(lineList)
        master_reactions=ct.Solution('masterTemp.cti')
        master_rxn_eqs=master_reactions.reaction_equations
    original_rxn_eqs=[]
    for i in np.arange(original_mechanism.n_reactions):
        if master_index[i]:
            NewModel.add_reaction(original_mechanism.reaction(i))
            original_rxn_count+=1
            original_rxn_eqs.append(original_mechanism.reaction_equation(i))
    if master_rxns:
        for i in np.arange(master_reactions.n_reactions):
            NewModel.add_reaction(master_reactions.reaction(i))
        
    if x=={}:
        for j in np.arange(original_rxn_count):
           
           if master_index[j]:
               if 'ThreeBodyReaction' in str(type(original_mechanism.reaction(j))):
                   NewModel.reaction(j).rate=original_mechanism.reaction(j).rate
               elif 'ElementaryReaction' in str(type(original_mechanism.reaction(j))):
                   NewModel.reaction(j).rate=original_mechanism.reaction(j).rate
               elif 'FalloffReaction' in str(type(original_mechanism.reaction(j))):
                   NewModel.reaction(j).high_rate=original_mechanism.reaction(j).high_rate
                   NewModel.reaction(j).low_rate=original_mechanism.reaction(j).low_rate
                   
                   if original_mechanism.reaction(j).falloff.type=='Troe':
                       NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                   if original_mechanism.reaction(j).falloff.type=='Sri':
                       NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
               elif 'ChemicallyActivatedReaction' in str(type(original_mechanism.reaction(j))):
                   NewModel.reaction(j).high_rate=original_mechanism.reaction(j).high_rate
                   NewModel.reaction(j).low_rate=original_mechanism.reaction(j).low_rate
                   if original_mechanism.reaction(j).falloff.type=='Troe':
                       NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                   if original_mechanism.reaction(j).falloff.type=='Sri':
                       NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
               elif 'PlogReaction' in str(type(original_mechanism.reaction(j))):
                   NewModel.reaction(j).rates=original_mechanism.reaction(j).rates
               elif 'ChebyshevReaction' in str(type(original_mechanism.reaction(j))):
                   NewModel.reaction(j).set_parameters(original_mechanism.reaction(j).Tmin,original_mechanism.reaction(j).Tmax,original_mechanism.reaction(j).Pmin,original_mechanism.reaction(j).Pmax,original_mechanism.reaction(j).coeffs)
    if x!={}:
        for j in np.arange(original_mechanism.n_reactions):
           if master_index[j]:
               try:
                   if 'ThreeBodyReaction' in str(type(original_mechanism.reaction(j))):
                       A=original_mechanism.reaction(j).rate.pre_exponential_factor
                       n=original_mechanism.reaction(j).rate.temperature_exponent
                       Ea=original_mechanism.reaction(j).rate.activation_energy
                       NewModel.reaction(j).rate=ct.Arrhenius(A*np.exp(x['r'+str(j)]['A']),n+x['r'+str(j)]['n'],Ea+x['r'+str(j)]['Ea'])
                   elif 'ElementaryReaction' in str(type(original_mechanism.reaction(j))):
                       A=original_mechanism.reaction(j).rate.pre_exponential_factor
                       n=original_mechanism.reaction(j).rate.temperature_exponent
                       Ea=original_mechanism.reaction(j).rate.activation_energy
                       NewModel.reaction(j).rate=ct.Arrhenius(A*np.exp(x['r'+str(j)]['A']),n+x['r'+str(j)]['n'],Ea+x['r'+str(j)]['Ea'])
                   elif 'FalloffReaction' in str(type(original_mechanism.reaction(j))):
                       A=original_mechanism.reaction(j).high_rate.pre_exponential_factor
                       n=original_mechanism.reaction(j).high_rate.temperature_exponent
                       Ea=original_mechanism.reaction(j).high_rate.activation_energy
                       NewModel.reaction(j).high_rate=ct.Arrhenius(A*np.exp(x['r'+str(j)]['A']),n+x['r'+str(j)]['n'],Ea+x['r'+str(j)]['Ea'])
                       A=original_mechanism.reaction(j).low_rate.pre_exponential_factor
                       n=original_mechanism.reaction(j).low_rate.temperature_exponent
                       Ea=original_mechanism.reaction(j).low_rate.activation_energy
                       NewModel.reaction(j).low_rate=ct.Arrhenius(A*np.exp(x['r'+str(j)]['A']),n+x['r'+str(j)]['n'],Ea+x['r'+str(j)]['Ea'])
                       if original_mechanism.reaction(j).falloff.type=='Troe':
                           NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                       if original_mechanism.reaction(j).falloff.type=='Sri':
                           NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                   elif 'ChemicallyActivatedReaction' in str(type(original_mechanism.reaction(j))):
                       A=original_mechanism.reaction(j).high_rate.pre_exponential_factor
                       n=original_mechanism.reaction(j).high_rate.temperature_exponent
                       Ea=original_mechanism.reaction(j).high_rate.activation_energy
                       NewModel.reaction(j).high_rate=ct.Arrhenius(A*np.exp(x['r'+str(j)]['A']),n+x['r'+str(j)]['n'],Ea+x['r'+str(j)]['Ea'])
                       A=original_mechanism.reaction(j).low_rate.pre_exponential_factor
                       n=original_mechanism.reaction(j).low_rate.temperature_exponent
                       Ea=original_mechanism.reaction(j).low_rate.activation_energy
                       NewModel.reaction(j).low_rate=ct.Arrhenius(A*np.exp(x['r'+str(j)]['A']),n+x['r'+str(j)]['n'],Ea+x['r'+str(j)]['Ea'])
                       if original_mechanism.reaction(j).falloff.type=='Troe':
                           NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                       if original_mechanism.reaction(j).falloff.type=='Sri':
                           NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                   elif 'PlogReaction' in str(type(original_mechanism.reaction(j))):
                       NewModel.reaction(j).rates=original_mechanism.reaction(j).rates
                   elif 'ChebyshevReaction' in str(type(original_mechanism.reaction(j))):
                       NewModel.reaction(j).set_parameters(original_mechanism.reaction(j).Tmin,original_mechanism.reaction(j).Tmax,original_mechanism.reaction(j).Pmin,original_mechanism.reaction(j).Pmax,original_mechanism.reaction(j).coeffs)
               except:
                   if 'ThreeBodyReaction' in str(type(original_mechanism.reaction(j))):
                       NewModel.reaction(j).rate=original_mechanism.reaction(j).rate
                   elif 'ElementaryReaction' in str(type(original_mechanism.reaction(j))):
                       NewModel.reaction(j).rate=original_mechanism.reaction(j).rate
                   elif 'FalloffReaction' in str(type(original_mechanism.reaction(j))):
                       NewModel.reaction(j).high_rate=original_mechanism.reaction(j).high_rate
                       NewModel.reaction(j).low_rate=original_mechanism.reaction(j).low_rate
                       if original_mechanism.reaction(j).falloff.type=='Troe':
                           NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                       if original_mechanism.reaction(j).falloff.type=='Sri':
                           NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                   elif 'ChemicallyActivatedReaction' in str(type(original_mechanism.reaction(j))):
                      NewModel.reaction(j).high_rate=original_mechanism.reaction(j).high_rate
                      NewModel.reaction(j).low_rate=original_mechanism.reaction(j).low_rate
                      if original_mechanism.reaction(j).falloff.type=='Troe':
                          NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                      if original_mechanism.reaction(j).falloff.type=='Sri':
                          NewModel.reaction(j).falloff=original_mechanism.reaction(j).falloff
                   elif 'PlogReaction' in str(type(original_mechanism.reaction(j))):
                      NewModel.reaction(j).rates=original_mechanism.reaction(j).rates
                   elif 'ChebyshevReaction' in str(type(original_mechanism.reaction(j))):
                      NewModel.reaction(j).set_parameters(original_mechanism.reaction(j).Tmin,original_mechanism.reaction(j).Tmax,original_mechanism.reaction(j).Pmin,original_mechanism.reaction(j).Pmax,original_mechanism.reaction(j).coeffs)
                   
               
    
    new_file=ctiw.write(NewModel)
    return new_file,original_rxn_eqs,master_rxn_eqs