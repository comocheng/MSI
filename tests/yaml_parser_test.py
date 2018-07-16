import sys
sys.path.append('.') #get rid of this at some point with central test script or when package is built

import MSI.simulations.instruments.shock_tube as st
import MSI.cti_core.cti_processor as pr
import MSI.simulations.yaml_parser as yp
import cantera as ct

parser = yt.Parser()
config = parser.load_to_obj('MSI/data/test_data/Hong_5.yaml')
loaded_tube = parser.parse_shock_tube_obj(loaded_exp=config)
print(loaded_tube)
