import cantera as ct
from . import soln2cti
#class holding different methods of processing cti files

class Processor(object): #handles one optimization but may add support for multiple later

    #only thing all processing needs is the cti file
    def __init__(self, path):
        self.cti_path = path
        self.solution = ct.Solution(path) #cantera handles its own errors, no need to except
        self.active_parameter_dictionary = {} # nothing for right now, options later for adding
    
    def add_active_parameter(self, r_index=-1,r_type='',dels=[], h_dels=[], l_dels=[],rate_list=[]):

        #checks so can't screw up input ie not give needed dels for a particular reaction type
        #only checks agains Reaction,ElementaryReaction, and Falloff, Plog will check in future

        if r_type=='EmentaryReaction':
            print('nothing')

        self.active_parameter_dictionary[r_index-1]=active_parameter(r_type,
                                                                     dels,h_dels,l_dels,
                                                                     rate_list)
    #sets default parameters for all reactions in solution according to corropsonding r type
    def set_default_parameters(self):
        for i in range(0, self.solution.n_reactions):
            r_type = self.soluion.reaction_type(i)
            if r_type == 'Reaction':
                print('x')
            elif r_type == 'ElementaryReaction':
                print('y')
            elif r_type == 'FalloffReaction':
                print('z')
            elif r_type == 'PlogReaction':
                print('a')

    def write_to_file(self,new_path=''):
        if new_path == '':
            new_path=self.cti_path.split(".cti")[0]+"_processed.cti"
            soln2cti.write(self.solution, new_path)
        else:
            soln2cti.write(self.solution,new_path)
        
        self.cti_path=new_path 
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
        self.hdels      = h_dels
        self.l_dels     = l_dels
        self.rate_list  = rate_list

