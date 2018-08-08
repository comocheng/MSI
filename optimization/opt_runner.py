import MSI.simulations as sim
import MSI.cti_core.cti_processor as pr
import MSI.optimization.matrix_loader as ml
import MSI.simulations.absorbance.curve_superimpose as csp
import MSI.simulations.yaml_parser as yp


# takes the data from one experiment and puts it in a dict of dicts
# that follows the format of the S matrix
# only need the last 3 elements of the interpolated absorbance
# absorbance is of form [interp_original,interp_abs_kinetic_sens,interp_abs_phys_sens,interp_abs_coef_sens]
# where each list element is a dict. keys are wavelengths, values are the sensitivities for that wavelength 
# psens should match interpolated_tp and species sens in size, but again is dict with wavelength keys
# index from 1, so if you have 3 experiments, their indices will be 1,2,3
def build_single_exp_dict(exp_index:int,
                          simulation:sim.instruments.shock_tube.shockTube,
                          interpolated_kinetic_sens:dict,
                          interpolated_tp_sens:list,
                          interpolated_species_sens:list,
                          interpolated_absorbance:list=[]):
    exp_dict = {}
    exp_dict['index']              = exp_index
    exp_dict['simulation']         = simulation
    exp_dict['ksens']              = interpolated_kinetic_sens
    exp_dict['temperature']        = interpolated_tp_sens[0]
    exp_dict['pressure']           = interpolated_tp_sens[1]
    exp_dict['species']            = interpolated_species_sens
    exp_dict['observables']        = simulation.observables
    if len(interpolated_absorbance) != 0:
        exp_dict['absorbance_ksens']   = interpolated_absorbance[1]
        exp_dict['absorbance_psens']   = interpolated_absorbance[2]
        exp_dict['perturbed_coef']     = interpolated_absorbance[3]
    return exp_dict


def load_exp_from_file(self,yaml_exp_file_list = []):
    for file in yaml_exp_file_list:
        continue
