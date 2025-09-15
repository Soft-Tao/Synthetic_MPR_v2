import numpy as np
import matplotlib.pyplot as plt
from components import Target, Aperture, Magnets, Focalplane
from processes import Beam
import ENDF_differential_cross_section
import contextlib
import os

class MPR:
    '''
    An integrated MPR system initialized by specifying the 4 components.
    It has the following functions:
    1. calculate l-E relation, energy resolution in a given energy range.
    2. find optimal straight focal plane position (position and tilt angle).
    3. calculate response matrix of the system.
    '''
    def __init__(self, target: Target, aperture: Aperture, magnet: Magnets, focalplane = None):
        self.target = target
        self.aperture = aperture
        self.magnet = magnet
        self.focalplane = focalplane
        self.resolution = None

    def performance(self, E_range = np.linspace(9, 21, 13), N_part = 50000, plot_show = True, plot_save = False):
        if self.focalplane is None:
            raise Exception("focalplane is not specified! Please specify focalplane first or run optimal_focalplane(). to find one.")
        output = []
        cross_section = ENDF_differential_cross_section.differential_cross_section(
            'E4R84432_e4.endf.endf2gnd.endf')

        for i, ene in enumerate(E_range):
            print(f"[{i}/{len(E_range)}] energy: {ene:.1f} MeV")
            beam = Beam.generate(self.target, self.aperture, ene, N_part, cross_section=cross_section)
            beam_transported = beam.trans(self.magnet)
            record = beam_transported.hit(self.focalplane)
            output.append([ene, record.l_mean, record.std_dev_l])
        
        output = np.array(output)
        if plot_show:
            plt.figure()
            plt.plot(output[:, 0], output[:, 1], '-s', color='k')
            plt.scatter(output[:, 0], output[:, 1], s=2, color='k')
            plt.title('L-E relation')
            plt.xlabel('energy [MeV]')
            plt.ylabel(r'$\bar{l}$ [m]')
            if plot_save:
                plt.savefig('L_E_relation.png', dpi=300, bbox_inches='tight')
            plt.show()
        
        resolution = output[1:-1, 2]/(output[2:, 1] - output[:-2, 1])*(output[2:, 0] - output[:-2, 0])/E_range[1:-1]*100
        if plot_show:
            plt.figure()
            plt.plot(E_range[1:-1], resolution, '-s', color='k')
            plt.scatter(E_range[1:-1], resolution, s=2, color='k')
            plt.title('Energy resolution') 
            plt.xlabel('energy [MeV]')
            plt.ylabel(r'Energy resolution [%]')
            plt.ylim(ymin=0)
            if plot_save:
                plt.savefig('energy_resolution.png', dpi=300, bbox_inches='tight')
            plt.show()
        
        self.resolution = resolution

    def optimal_focalplane(self, target: Target, aperture: Aperture, magnet: Magnets, merit_weight: str, fp_position: list, fp_angle: list, merit_weight_lst = None):
        '''
        About the weight:
        Its a zero-dimensional number to represent an average energy resolution of the system.
        If merit_weight = 'uniform', the weight of each energy is 1.
        Otherwise, you have to specify a merit_weight_lst (len=11) to Energy_lst (np.linspace(10, 20, 11)).
        '''
        print("Starting to find optimal focal plane position..." \
        "Notice: self.focalplane will be replaced by optimal focal plane.")
        if merit_weight == 'uniform':
            merit_weight_lst = np.ones(11) # uniform weight, default length of resolution list is 11
        else:
            if merit_weight_lst is None:
                raise Exception("merit_weight_lst is required when merit_weight is not [uniform]!!!") 

        merit_matrix = np.zeros((len(fp_position), len(fp_angle)))
        for i, position in enumerate(fp_position):
            for j, angle in enumerate(fp_angle):
                print(f"\rProcess: {100*(i*len(fp_angle)+j)/(len(fp_position)*len(fp_angle)):.2f}%", end='')
                geometry = [(x, position + x * np.tan(angle)) for x in np.linspace(-0.5, 0.5, 11)]
                self.focalplane = Focalplane('arbitrary', position=position, geometry=geometry)

                with contextlib.redirect_stdout(open(os.devnull, 'w')):
                    self.performance(N_part=5000, plot_show=False)
                merit_matrix[i, j] = np.sum(self.resolution[:][1]*merit_weight_lst)/np.sum(merit_weight_lst)
        
        min_index = np.unravel_index(np.argmin(merit_matrix), merit_matrix.shape)
        print('\n')
        print(f"Optimization of the focal plane is completed. Minimum merit value is {merit_matrix[min_index[0]][min_index[1]]:.4f}." \
        f"The optimal focal plane position is {fp_position[min_index[0]]: .4f}, tilt angle with 'normal' direction is {fp_angle[min_index[1]]: .4f}.")
        geometry_optimal = [(x, fp_position[min_index[0]] + x * np.tan(fp_angle[min_index[1]])) for x in np.linspace(-0.5, 0.5, 11)]
        self.focalplane = Focalplane('arbitrary', position=fp_position[min_index[0]], geometry=geometry_optimal)

target = Target(type='CH2', mass_thickness=15, shape="circle", area=0.001, H_W_ratio=1.0)

aperture = Aperture(distance=0.2, shape="circle", solid_angle=40, H_W_ratio=1.0)

Magnet = Magnets('beam_trans_benchmark/TM.txt', reference_energy=14)

mpr = MPR(target, aperture, Magnet)

mpr.optimal_focalplane(target, aperture, Magnet, merit_weight='uniform', fp_position=np.linspace(0, 0.5, 8), fp_angle=np.linspace(10, 60, 25))
mpr.performance(plot_show=True, N_part=5000)