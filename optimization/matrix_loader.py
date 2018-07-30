import numpy as np
import pandas as pd

class OptMatrix(object):
    def __init__(self):
        self.matrix = None
 
    #loads one experiment into self.matrix. Decides padding based on previous matrix or handle based on total exp num?
    def load_S(self, exp_num, list_of_interpolated_kinetic_sens,
                              list_of_interpolated_tp_sens,
                              list_of_interpolated_species_sens,
                              list_of_interpolated_absorbance_sens,
                              observables_list):
    
        for i in range(0,exp_num):
            full_exp_matrix = self.load_single_exp(i,list_of_interpolated_kinetic_sens[i],
                                                     list_of_interpolated_tp_sens[i],
                                                     list_of_interpolated_species_sens[i],
                                                     list_of_interpolated_absorbance_sens[i],
                                                     observables_list)
    
    def load_single_exp(self, exp_num, interpolated_kinetic_sens,
                                       interpolated_tp_sens,
                                       interpolated_species_sens, 
                                       observables_list):
        '''
        interpolated_kinetic_sens: list of 3 numpy arrays. Each array contains sheets representing one observable
                                   interpolated from experimental data
        interpolated tp sens: interpolated temperature and pressure sensitivities from some experimental data
        interpolated_species_sens: species sensitiviites listed by observables
        '''
        '''
        obs_matrix = None
        base_matrix = interpolated_kinetic_sens[0][0] #make a base to build the matrix off of
        for i,obs in enumerate(observables_list):
            #kinetic sensitiviities
            #build a long horizontal array then stack
            print(interpolated_kinetic_sens[0][i].shape)
            base_matrix=np.hstack((interpolated_kinetic_sens[0][i], 
                                interpolated_kinetic_sens[1][i],
                                interpolated_kinetic_sens[2][i],))
                                #interpolated_tp_sens[0][obs].values.reshape(50,1),
                                #interpolated_tp_sens[1][obs].values.reshape(50,1)))
            #for spec_sen in interpolated_species_sens:
            #    base_matrix=np.hstack((base_matrix,spec_sen))
            # wrong way to loop
            if i == 0:
                obs_matrix = base_matrix
            else:
                obs_matrix = np.vstack((self.matrix,base_matrix)) 
        #need to do abs, and handle padding,
        #do calc to figure out size of the 0 chunks for mult exp
        return self.matrix
        '''
        pass
    def load_Y():
        pass

    def build_Z():
        pass
