import cantera as ct
import soln2cti
import sys
#class holding different methods of processing cti files

class Processor(object):

    #only thing all processing needs is the cti file
    def __init__(self, path):
        self.cti_path = path
        self.solution = ct.Solution(path) #cantera handles its own errors, no need to except
    
    #varialble number of args, can take a file path, list of numbers or both
    #removes the specified reactions from the solution object created in constructor
    def prune(*args):
        if len(args) == 1:
            if isinstance(args[0],str):
                print("temp")
            elif isinstance(args[0],list):
                print("temp")

            else:
                print("When using a single argument, give only a file path or list of integers."
               
        elif len(args) == 2:
            if not isinstance(args[0],str) or not isinstance(args[1],list):
                print("Please enter parameters as prune(path,list) when using 2 arguments")
                return 

        else:
            print("Incorrect number of arguments.")
            return
