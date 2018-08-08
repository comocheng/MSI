import numpy as np
import pandas as pd

class OptMatrix(object):
    def __init__(self):
        self.matrix = None
 
    #loads one experiment into self.matrix. Decides padding based on previous matrix or handle based on total exp num?
    def load_S(self, exp_dict_list:list):
        #preprocessing for padding
        num_exp = len(exp_dict_list)
        pert_coef = {} #build a dict matching pert_coef to their experiment and wavelength.
                       #length of the dict gives padding information
        for exp in exp_dict_list:
            if 'perturbed_coef' not in exp.keys():
                continue
            perturbed_for_exp = exp['perturbed_coef']
            for x in perturbed_for_exp:
                if x[0][2] not in pert_coef.keys():
                    pert_coef[x[0][2]] = [x[1]]
                else:
                    pert_coef[x[0][2]].append(x[1])
        num_ind_pert_coef = len(pert_coef)
        #print(num_ind_pert_coef," sigmas")
        #establish # of independent pert before hand, to proper pad the observables, put in list, make a dict of cc,
        # values will be a list of tabs data?
        # use the list to get the padding size
        for i in range(0,num_exp):
            #for each observable
            exp_matrix = None
            for j,obs in enumerate(exp_dict_list[i]['observables']):
                #kinetic sensitiviities
                #build a long horizontal array then stack
                #print(exp_dict_list[i]['ksens'][0][i].shape)
                obs_matrix=np.hstack((exp_dict_list[i]['ksens']['A'][i], 
                                    exp_dict_list[i]['ksens']['N'][i],
                                    exp_dict_list[i]['ksens']['Ea'][i],))
                #for spec_sen in interpolated_species_sens:
                #    obs_matrix=np.hstack((obs_matrix,spec_sen))
                # wrong way to loop
                                
                #do the prepadding psens
                num_zeros = 0
                for jj in range(0,i): #loop the exp before exp i
                    #tp padding
                    num_zeros+=2 
                    #species padding
                    num_zeros+= len(exp_dict_list[jj]['observables'])
                obs_matrix = np.hstack((obs_matrix,
                                        np.zeros((obs_matrix.shape[0],num_zeros))))
                #do the psens
                print("ksens shape:", obs_matrix.shape,"\ntemp and pressure shape:",exp_dict_list[i]['temperature'][obs].shape)
                t_stack = np.vstack((exp_dict_list[i]['temperature'][obs],
                                    exp_dict_list[i]['pressure'][obs])).T
                print("stacked tp shape:",t_stack.shape,"\n")
                print(obs_matrix[:,0])
                print("what the fuck")
                print(t_stack[:,0])
                obs_matrix = np.hstack((obs_matrix,t_stack))
                                        #exp_dict_list[i]['temperature'][obs],
                                        #exp_dict_list[i]['pressure'][obs]))
                #psens post padding
                num_zeros = 0
                for jj in range(i+1,num_exp): 
                    #tp padding
                    num_zeros+=2 
                    #species padding
                    num_zeros+= len(exp_dict_list[jj]['observables'])
                obs_matrix = np.hstack((obs_matrix,
                                        np.zeros((obs_matrix.shape[0],num_zeros))))

                #sigma padding, how do
                
                #do the vertical stacking
                if i == 0:
                    exp_matrix = obs_matrix
                else:
                    exp_matrix = np.vstack((self.matrix,obs_matrix)) 
                
            #loop the abs data
                #do abs ksens


                #pre padding

                #psens

                #postpadding

                #sigmas and their padding

                #vstack

            #vstack the whole exp

        #vstack the indentity

        #maybe have metadata in the matrix or in some other dict, to give info on where and how get data out?
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
        obs_matrix = interpolated_kinetic_sens[0][0] #make a base to build the matrix off of
        for i,obs in enumerate(observables_list):
            #kinetic sensitiviities
            #build a long horizontal array then stack
            print(interpolated_kinetic_sens[0][i].shape)
            obs_matrix=np.hstack((interpolated_kinetic_sens[0][i], 
                                interpolated_kinetic_sens[1][i],
                                interpolated_kinetic_sens[2][i],))
                                #interpolated_tp_sens[0][obs].values.reshape(50,1),
                                #interpolated_tp_sens[1][obs].values.reshape(50,1)))
            #for spec_sen in interpolated_species_sens:
            #    obs_matrix=np.hstack((obs_matrix,spec_sen))
            # wrong way to loop
            if i == 0:
                obs_matrix = obs_matrix
            else:
                obs_matrix = np.vstack((self.matrix,obs_matrix)) 
        #need to do abs, and handle padding,
        #do calc to figure out size of the 0 chunks for mult exp
        return self.matrix
        '''
        pass
    def load_Y():
        pass

    def build_Z():
        pass
