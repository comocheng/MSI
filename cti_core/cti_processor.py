import cantera as ct
from . import soln2cti
#class holding different methods of processing cti files

class Processor(object): #handles one optimization but may add support for multiple later

    #only thing all processing needs is the cti file
    def __init__(self, path, p_flag = 0):
        self.cti_path = path
        self.solution = ct.Solution(path) #cantera handles its own errors, no need to except
        self.active_parameter_dictionary = {} # nothing for right now, options later for adding
        self.valid_reactions = ['ElementaryReaction',
                                'ThreeBodyReaction',
                                'FalloffReaction',
                                'PlogReaction']
        self.param_path = ''#may give option to load in a param file
        if p_flag == 1:
            self.set_default_parameters()
    #v flag for verbose output, may add out for A=, n= , Ea = etc later, and rate list for plog support
    def get_active_parameter(self,r_index, verbose=0):
        if len(self.active_parameter_dictionary) == 0:
            print("Error: active parameter dictionary empty. Please initialize or add a parameter to the dictionary")
            return -1
        elif r_index-1 < 0:
            print("Error: negative reaction index")
            return -1
        elif r_index-1 >= self.solution.n_reactions:
            print("Error: invalid index out of bounds")
            return -1
        elif r_index not in self.active_parameter_dictionary.keys():
            print("Error: reaction {0} at index{1} has no parameter entry".format(self.solution.reaction(i-1),i))
            return -1
        if v==1:
            print("Reaction {0}:\nType: {1}\ndels: {2}\nh_dels: {3}\nl_dels: {4}\n rate_list: {5}".format(
                self.solution.reaction(r_index - 1),
                self.active_parameter_dictionary[r_index-1].r_type,
                self.active_parameter_dictionary[r_index-1].dels,
                self.active_parameter_dictionary[r_index-1].h_dels,
                self.active_parameter_dictionary[r_index-1].l_dels,
                self.active_parameter_dictionary[r_index-1].rate_list))
        return self.active_parameter_dictionary[r_index-1]
    
    #processor maintains the activate parameters, they are not manipulated in this class
    def add_active_parameter(self, r_index=-1,r_type='',dels=[], h_dels=[], l_dels=[],rate_list=[]):

        #checks so can't screw up input ie not give needed dels for a particular reaction type
        #only checks agains Reaction,ElementaryReaction, and Falloff, Plog will check in future
        
        if r_type not in self.valid_reactions:
            print('Error: Non supported reaction type {0}'.format(r_type))
            print('Valid Reaction Types:', self.valid_reactions)
        if r_type=='ElementaryReaction' or r_type == 'ThreeBodyReaction':
            if len(h_dels)!=0 or len(l_dels)!=0 or len(rate_list)!=0:
                print('Error: Invalid parameter given, {0}  only takes the dels list'.format(r_type))
                return False
            if len(dels) == 0:
                print('Error: dels cannot be empty for generating ThreeBodyReaction parameter')
                return False
            if len(dels) != 3:
                print('Error: dels takes 3 arguments: A,n,Ea')
                return False
            for x in dels:
                if not isinstance(x,float):
                    print('Error: {0} is not a float, all del vals must be floats'.format(x))
                    return False
        if r_type=='FalloffReaction':
            if len(dels)!=0 or len(rate_list)!=0:
                print('Error: Invalid parameter given, {0}  only takes the h_dels and l_dels '.format(r_type))
                return False
            if len(h_dels) == 0:
                print('Error: h_dels cannot be empty for generating FalloffReaction parameter')
                return False
            if len(l_dels) == 0:
                print('Error: l_dels cannot be empty for generating FalloffReaction parameter')
                return False
            if len(h_dels) != 3:
                print('Error: h_dels takes 3 arguments: A,n,Ea')
                return False
            if len(l_dels) != 3:
                print('Error: l_dels takes 3 arguments: A,n,Ea')
                return False

            for x,y in zip(h_dels,l_dels):
                if not isinstance(x,float):
                    print('Error: {0} in h_dels is not a float, all del vals must be floats'.format(x))
                if not isinstance(y,float):
                    print('Error: {0} in l_dels is not a float, all del vals must be floats'.format(y))
                    return False

        elif r_type == 'PlogReaction':
            print('not supported yet')
            return False
        
        
        self.active_parameter_dictionary[r_index-1]=active_parameter(r_type,
                                                                     dels,h_dels,l_dels,
                                                                     rate_list)
        return True 
    #sets default parameters for all reactions in solution according to corropsonding r type
    def set_default_parameters(self):
        for i in range(0, self.solution.n_reactions):
            if isinstance(self.solution.reaction(i),ct._cantera.ElementaryReaction):
                self.add_active_parameter(r_index = i,r_type = 'ElementaryReaction',dels=[0.0,0.0,0.0])
            elif isinstance(self.solution.reaction(i),ct._cantera.ThreeBodyReaction):
                self.add_active_parameter(r_index = i,r_type = 'ThreeBodyReaction',dels=[0.0,0.0,0.0])
            elif isinstance(self.solution.reaction(i),ct._cantera.FalloffReaction):
                self.add_active_parameter(r_index = i,r_type = 'FalloffReaction',h_dels=[0.0,0.0,0.0],l_dels=[0.0,0.0,0.0])
            elif isinstance(self.solution.reaction(i),ct._cantera.PlogReaction):
                print('PlogReaction not supported yet, skipping')
            else:
                print('Unsupported Reaction Type {0},index {1} skipping'.format(self.solution.reaction(i).reaction_type(),i+1))
    #write the new cti file, note original file is always preserved unless its path is given as new_path
    #also note that _processed.cti will repeatably be rewritten if no new path is specified
    def write_to_file(self,new_path=''):
        if new_path == '':
            new_path=self.cti_path.split(".cti")[0]+"_processed.cti"
            soln2cti.write(self.solution, new_path)
        else:
            soln2cti.write(self.solution,new_path)
        
        self.cti_path=new_path 
        return new_path
    #write the active parameter information to file
    #naming scheme behaves in same way as write_to_file, may add option to do spefiic reactions, not only all
    def write_parameters_to_file(self,new_path=''):
        if new_path == '':
            new_path=self.cti_path.split(".cti")[0]+"_processed.param"
        f = open(new_path,'w')
        for r_index in self.active_parameter_dictionary.keys():
            data = "Reaction {0}:\nType: {1}\ndels: {2}\nh_dels: {3}\nl_dels: {4}\nrate_list: {5}\n".format(
            self.solution.reaction(r_index + 1),
            self.active_parameter_dictionary[r_index].r_type,
            self.active_parameter_dictionary[r_index].dels,
            self.active_parameter_dictionary[r_index].h_dels,
            self.active_parameter_dictionary[r_index].l_dels,
            self.active_parameter_dictionary[r_index].rate_list)
            f.write(data)
        self.param_path=new_path
        return new_path

    #expects a list of integers representing reaction indices to remove
    #Cantera as of 2.3 does not have a native remove reaction function
    #assumes input indices are from 1 to n, then subs for cantera to do 0 to n-1
    def remove_reactions(self, to_remove:list):
        clean_reactions=[]
        for i in to_remove:
            if not isinstance(i,int):
                print("{0} not an integer, will not be removed".format(i))
            else:
                print("remove index {0}, reaction {1}".format(i,self.solution.reaction(i-1)))

        for i in range(0,self.solution.n_reactions):
            if i not in to_remove:
                clean_reactions.append(self.solution.reaction(i))
                  
        self.solution = ct.Solution(thermo='IdealGas',
                                    kinetics='GasKinetics',
                                    species=self.solution.species(),
                                    reactions=clean_reactions)

    #appends list of reaction indices with those from a file
    #expects format of one index per line, more behavior in future?
    def append_list(self, path, list_to_add):
       try:
           f = open(path,'r')
           for i,line in enumerate(f):
               try:
                  list_to_add.append(int(line.strip()))

               except ValueError as e:
                  print('Error on index {0}: {1}\n Skipping index'.format(line,e))
                
       except IOError as e:
            print("Error: {0}".format(e))
       return list_to_add

    #variable number of args, can take a file path, list of numbers or both, or single index
    #removes the specified reactions from the solution object created in constructor
    def prune(self,*args):
        if len(args)==1:
            if isinstance(args[0],str):
                l = self.append_list(args[0], [])
                self.remove_reactions(l)
            elif isinstance(args[0],list):
                self.remove_reactions(args[0])
            elif isinstance(args[0], int):
                self.remove_reactions([args[0]])
            else:
                print("When using a single argument, give only a file path", 
                       "or list of integers, or a single integer")
        elif len(args)==2:
             if not isinstance(args[0],str) or not isinstance(args[1],list):
                print("Please enter parameters as prune(path,list) when using 2 arguments")
                return 
             else:
                l = self.append_list(args[0],args[1])
                self.remove_reactions(l)
        else:
            print("Incorrect number of arguments.")

#represents an active parameter set which applies to a set of reactions
class active_parameter(object):
    def __init__(self, r_type='',dels=[], h_dels=[], l_dels=[],rate_list=[]):
        self.r_type     = r_type
        self.dels       = dels
        self.h_dels      = h_dels
        self.l_dels     = l_dels
        self.rate_list  = rate_list

