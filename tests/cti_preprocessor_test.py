import importlib.util
spec = importlib.util.spec_from_file_location("cti_preprocessor", "../cti_preprocessing/soln2cti.py")
cti_pp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cti_pp)

import cantera as ct
gas = ct.Solution("../data/test_data/lam.cti")

cti_pp.write(gas, "../data/test_data/cti_test.cti")
