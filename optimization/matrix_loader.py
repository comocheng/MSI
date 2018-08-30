import numpy as np
import pandas as pd

class OptMatrix(object):
    def __init__(self):
        self.S = None
 
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
                if j < len(exp_dict_list[i]['ksens']['A']):
                    obs_matrix=np.hstack((exp_dict_list[i]['ksens']['A'][j], 
                                        exp_dict_list[i]['ksens']['N'][j],
                                        exp_dict_list[i]['ksens']['Ea'][j],))
                                    
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
                    t_stack = np.vstack((exp_dict_list[i]['temperature'][obs],
                                        exp_dict_list[i]['pressure'][obs])).T
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
                    obs_matrix = np.hstack((obs_matrix,np.zeros((obs_matrix.shape[0],num_ind_pert_coef))))                    
                    #do the vertical stacking
                    if j == 0:
                        exp_matrix = obs_matrix
                    else:
                        exp_matrix = np.vstack((exp_matrix,obs_matrix)) 
                else:
                    break
            #loop the abs data, have to do it by wavelength, build the sigmas at once sep?
                #do abs ksens


                #pre padding

                #psens

                #postpadding

                #sigmas and their padding

                #vstack

            #vstack the whole exp

        #vstack the indentity

        #maybe have metadata in the matrix or in some other dict, to give info on where and how get data out?
    def load_Y(self, exp_dict_list:list):
        def natural_log_difference(experiment,model):
            natural_log_diff = np.log(experiment) - np.log(model)
            return natural_log_diff
        
        Y = []
        YdataFrame = []
        for i,exp_dic in enumerate(exp_dict_list):
            for j,observable in enumerate((exp_dic['mole_fraction_observables']+
                                           exp_dic['concentration_observables'])):
                natural_log_diff = natural_log_difference(exp_dic['experimental_data'][j][observable].values,
                                                          exp_dic['simulation'].timeHistoryInterpToExperiment[observable].dropna().values)
                natural_log_diff =  natural_log_diff.reshape((natural_log_diff.shape[0],
                                                  1))
                
                tempList = [observable+'_'+'experiment'+str(i)]*np.shape(natural_log_diff)[0]
                YdataFrame.extend(tempList)
                
                Y.append(natural_log_diff)
            if 'absorbance_observables' in list(exp_dic.keys()):
                wavelengths = list(exp_dic['absorbance_ksens'].keys())
                
                for k,wl in enumerate(wavelengths):
                    natural_log_diff = natural_log_difference(exp_dic['absorbance_experimental_data'][k]['Absorbance_'+str(wl)].values,exp_dic['absorbance_model_data'][wl].values)
                    natural_log_diff =  natural_log_diff.reshape((natural_log_diff.shape[0],
                                                  1))
                    
                    tempList = [str(wl)+'_'+'experiment'+str(i)]*np.shape(natural_log_diff)[0]
                    YdataFrame.extend(tempList)
                    
                    
                    Y.append(natural_log_diff)
        
        Y = np.vstack((Y))
        YdataFrame = pd.DataFrame({'value': YdataFrame,'ln_difference': Y})
        
        return Y, YdataFrame

    def build_Z():
        pass
