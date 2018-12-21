import numpy as np

class testing_code(object):
    def __init__(self,shape_of_X,value=.1,Y_matrix = None):
        self.S_matrix = None
        self.s_matrix = None
        self.Y_matrix = None
        self.y_matrix = None
        self.shape_of_X = shape_of_X
        self.value_to_perturb_by = value
        
        
        
        
        
    def perturb_X(self,counter):
        
       
        zeroed_out = np.zeros((self.shape_of_X))
        zeroed_out[counter,:] = self.value_to_perturb_by
        
        return zeroed_out
    
    
    
    def calculate_Sij(self,counter,lengths_of_experimental_data_list):
        start=0 
  
        S_residuals = np.zeros((np.shape(self.S_matrix)))
        for x in range(len(lengths_of_experimental_data_list)):
            for y in range(len(lengths_of_experimental_data_list[x])):
                stop = self.simulation_lengths_of_experimental_data[x][y] + start
                original = self.S_matrix[start:stop,counter]
                numerator = self.Y_matrix[start:stop,:]
                denominator = self.value_to_perturb_by
                Sij = np.divide(numerator,denominator)
                residuals = original-Sij
                S_residuals[start:stop,counter] = residuals
                start = start + self.simulation_lengths_of_experimental_data[x][y]
        
        return S_residuals
                
     