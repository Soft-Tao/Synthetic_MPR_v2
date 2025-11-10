import numpy as np
import matplotlib.pyplot as plt
from components import Target, Aperture, Magnets, Focalplane
from processes import Beam
import ENDF_differential_cross_section
import contextlib
import os
from tqdm import tqdm

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
            raise Exception("focalplane is not specified! Please specify focalplane first or run optimal_focalplane() to find one.")
        output = []
        efficiency_lst = []
        cross_section = ENDF_differential_cross_section.differential_cross_section(
            'E4R84432_e4.endf.endf2gnd.endf')
        if self.magnet.type == '2d':
            if plot_show: 
                plt.figure(figsize=(10, 4),dpi = 140)
                plt.xlabel("l [m]")
                plt.ylabel("y [m]")
            for i, ene in enumerate(E_range):
                print(f"[{i}/{len(E_range)}] energy: {ene:.1f} MeV")
                beam = Beam.generate(self.target, self.aperture, ene, N_part, cross_section=cross_section)
                beam_transported = beam.trans(self.magnet)
                efficiency_lst.append(beam.N_beampart / beam.N_neutrons)
                record = beam_transported.hit(self.focalplane)
                if plot_show: plt.scatter(record.l_hits, record.y_hits, s=0.1, label = f"{ene} MeV")
                output.append([ene, record.l_mean, record.std_dev_l])
            if plot_show: plt.legend(markerscale = 20)
            if plot_show: plt.show()
            
            output = np.array(output)
            efficiency_lst = np.array(efficiency_lst)
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
                fig, ax1 = plt.subplots()
                
                # Plot energy resolution on left y-axis
                color = 'k'
                ax1.set_xlabel('energy [MeV]')
                ax1.set_ylabel(r'Energy resolution [%]', color=color)
                ax1.plot(E_range[1:-1], resolution, '-s', color=color)
                ax1.scatter(E_range[1:-1], resolution, s=2, color=color)
                ax1.tick_params(axis='y', labelcolor=color)
                ax1.set_ylim(ymin=0)
                
                # Create second y-axis for efficiency
                ax2 = ax1.twinx()
                color = 'r'
                ax2.set_ylabel('Efficiency', color=color)
                ax2.plot(E_range[1:-1], efficiency_lst[1:-1], '-s', color=color)
                ax2.scatter(E_range[1:-1], efficiency_lst[1:-1], s=2, color=color)
                ax2.tick_params(axis='y', labelcolor=color)
                
                plt.title('Energy resolution & efficiency')
                if plot_save:
                    plt.savefig('energy_resolution.png', dpi=300, bbox_inches='tight')
                plt.show()
            
            self.resolution = resolution

    def optimal_focalplane(self, fp_position: list, fp_angle: list, merit_weight_lst = None, merit_weight: str = 'uniform'):
        '''
        About the weight:
        Merit is a zero-dimensional number to represent an average energy resolution of the system. Weight is what we need when calculating the average.
        If merit_weight = 'uniform', the weight of each energy is 1.
        Otherwise, you have to specify a merit_weight_lst (len=11) to Energy_lst (np.linspace(10, 20, 11)).
        =====2d=====
        fp_postion: list of z (in beam reference)
        =====2d=====
        fp_postion: list of (x, y) tuple (in 3d xyz coordinate system)
        '''
        print("Starting to find optimal focal plane position..." \
        "Notice: self.focalplane will be replaced by optimal focal plane.")
        if merit_weight == 'uniform':
            merit_weight_lst = np.ones(11) # uniform weight, default length of resolution list is 11
        else:
            if merit_weight_lst is None:
                raise Exception("merit_weight_lst is required when merit_weight is not [uniform]!!!") 
            else:
                merit_weight_lst = np.array(merit_weight_lst)

        # calculate E_range = np.linspace(9, 21, 13), N_part = 5000, store total beam_out @ magnet's exit.
        cross_section = ENDF_differential_cross_section.differential_cross_section(
            'E4R84432_e4.endf.endf2gnd.endf')
        if self.magnet.type == '2d':
            beam_lst = []
            for i,ene in enumerate(np.linspace(9, 21, 13)):
                print(f"[{i}/{len(np.linspace(9, 21, 13))}] energy: {ene:.1f} MeV")
                beam = Beam.generate(self.target, self.aperture, ene, 5000, cross_section=cross_section)
                beam_transported = beam.trans(self.magnet)
                beam_lst.append(beam_transported)
            
            merit_matrix = np.zeros((len(fp_position), len(fp_angle)))
            for i, position in enumerate(fp_position):
                for j, angle in enumerate(fp_angle):
                    print(f"\rProcess: {100*(i*len(fp_angle)+j)/(len(fp_position)*len(fp_angle)):.2f}%", end='')
                    geometry = [(x, position + x * np.tan(angle*np.pi/180)) for x in np.linspace(-0.5, 0.5, 11)]
                    self.focalplane = Focalplane('arbitrary', position=position, geometry=geometry)
                    output = []
                    for ene, beam_out in zip(np.linspace(9, 21, 13), beam_lst):
                        record = beam_out.hit(self.focalplane)
                        output.append([ene, record.l_mean, record.std_dev_l])
                    output = np.array(output)
                    self.resolution = output[1:-1, 2]/(output[2:, 1] - output[:-2, 1])*(output[2:, 0] - output[:-2, 0])/np.linspace(9, 21, 13)[1:-1]*100
                    merit_matrix[i, j] = np.sum(self.resolution[:][1]*merit_weight_lst)/np.sum(merit_weight_lst)
            
            min_index = np.unravel_index(np.argmin(merit_matrix), merit_matrix.shape)
            print('\n')
            print(f"Optimization of the focal plane is completed. Minimum merit value is {merit_matrix[min_index[0]][min_index[1]]:.4f}." \
            f"The optimal focal plane position is {fp_position[min_index[0]]: .4f}, tilt angle with 'normal' direction is {fp_angle[min_index[1]]: .4f}.")
            geometry_optimal = [(x, fp_position[min_index[0]] + x * np.tan(fp_angle[min_index[1]]*np.pi/180)) for x in np.linspace(-0.5, 0.5, 101)]
            self.focalplane = Focalplane('arbitrary', position=fp_position[min_index[0]], geometry=geometry_optimal)
    
    def save_focalplane(self, save_path: str):
        if self.focalplane is None:
            raise Exception("focalplane is not specified! Please specify focalplane first or run optimal_focalplane() to find one.")
        else:
            with open(os.path.join(save_path, "fp.txt"), 'w') as f:
                for line in self.focalplane.geometry:
                    f.write(f"{line[0]} {line[1]}\n")
            f.close()
            print(f"current focalplane has been successfully saved as: {os.path.join(save_path, 'fp.txt')}!")
    
    def response_matrix(self, E_min: float = 10, E_max: float = 20, sample_times:int = 100000 , N_part: int = 100, save_path = None, plot_save = False):
        '''
        Calculates the response_matrix of one MPR system. It's normalized to energy-depended efficiency.

        Proton strike position l: l_bins depends on the geometry of focalplane ((x, z) nodes will be regarded as bin edges).
        Parameters:
        E_min, E_max: minimun and maximun of the incoming nuetron energy.
        N_part: length of Beam object per sampling. (Total counts = sample_times * N_part) The final fluctuation of the response matrix would be decreased by increasing total_counts.
        '''
        cross_section = ENDF_differential_cross_section.differential_cross_section(
            'E4R84432_e4.endf.endf2gnd.endf')
        E_lst = []
        l_lst = []
        print("Start sampling...")
        for i in tqdm(range(sample_times)):
            E = np.random.uniform(E_min, E_max)
            E_lst += [E for x in range(N_part)]
            beam = Beam.generate(self.target, self.aperture, E, N_part, cross_section)
            beam_transported = beam.trans(self.magnet)
            record = beam_transported.hit(self.focalplane)
            l_lst += record.l_hits.tolist()
        print(f"Total counts of recorded protons: {sample_times * N_part}")
        plt.figure()
        H, l_edges, E_edges, _ = plt.hist2d(l_lst, E_lst, bins=(self.focalplane.geometry[:,2], np.linspace(E_min, E_max, 201)), cmap = 'jet')
        # normalize
        sum_col = np.sum(H, axis=0)
        sum_col[sum_col == 0] = 1
        H = H / sum_col
        # take efficiency into account
        efficiency_lst = []
        print("Calculating efficiency in each energy bin...")
        for i in tqdm(range(H.shape[1])):
            Ec = (E_edges[i] + E_edges[i+1]) / 2
            beam_ = Beam.generate(self.target, self.aperture, Ec, 50000, cross_section)
            efficiency_lst.append(beam_.N_beampart / beam_.N_neutrons)
        #efficiency_lst_reshaped = np.reshape(np.array(efficiency_lst), (-1, 1))
        H = H * np.array(efficiency_lst)
        if plot_save:
            if save_path is None:
                raise Exception("you CANNOT save the response matrix without specifying a save path!")
            else:
                np.savez(os.path.join(save_path, "response_matrix.npz"), H = H, l_edges = l_edges, E_edges = E_edges)
        aspect_ratio = (l_edges[-1] - l_edges[0]) / (E_edges[-1] - E_edges[0])

        ax = plt.imshow(H.T, origin='lower', extent=[l_edges[0], l_edges[-1], E_edges[0], E_edges[-1]], cmap='jet', aspect=aspect_ratio)
        plt.title('Response Matrix')
        plt.colorbar()
        plt.xlabel('Proton strike position [m]')
        plt.ylabel('Incoming neutron energy [MeV]')
        plt.show()
