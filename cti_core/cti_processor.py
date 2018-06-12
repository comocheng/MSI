import cantera as ct
import soln2cti
import sys
#class holding different methods of processing cti files

class Processor(object):

    #only thing all processing needs is the cti file
    def __init__(self, path):
        self.cti_path = path
        self.solution = ct.Solution(path) #cantera handles its own errors, no need to except
    
    def write_to_file(self):
        soln2cti(self.solution, self.path)
    
    #expects a list of integers representing reaction indices to remove
    def remove_reactions(self, to_remove):
        for i in to_remove:
            if not isinstance(i,int):
                print("invalid index type, skipping")
                continue
            elif i<0:
                print("invalid index, skipping")
                continue
            clean_reactions = []
            for i in range(0,self.solution.n_reactions):
                if i not in to_remove:
                    clean_reactions.append(self.solution.reaction(i))

            self.solution = ct.Solution(thermo='IdealGas',
                                        kinetics='GasKinetics',
                                        species=self.species(),
                                        reactions=clean_reactions)

    #appends list of reaction indices with those from a file
    #expects format of one index per line, more behavior in future?
    def append_list(self, path, list_to_add):
        try:
            f = open(path,'r')
            for i,line in enumerate(f):
                try:
                    list_to_add.append(int(line))
                except ValueError as e:
                    print("Error: {0}".format(e))
        except IOError as e:
            print("Error: {0}".format(e))
            return list_to_add

    #varialble number of args, can take a file path, list of numbers or both
    #removes the specified reactions from the solution object created in constructor
    def prune(self,*args):
        if len(args)==1:
            if isinstance(args[0],str):
                l = appendList(args[0], [])
                remove_reactions(l)
            elif isinstance(args[0],list):
                remove_reactions(list)
            else:
                print("When using a single argument, give only a file path", 
                       "or list of integers.")
        elif len(args)==2:
             if not isinstance(args[0],str) or not isinstance(args[1],list):
                print("Please enter parameters as prune(path,list) when using 2 arguments")
                return 
             else:
                l = appendList(args[0],args[1])
                remove_reactions(l)
        #only hits here if wrong number of arguments given
        print("Incorrect number of arguments.")
        return
