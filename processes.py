import copy, os
import matplotlib.pyplot as plt
import numpy as np
from components import Target, Aperture, Magnets, Focalplane
import ENDF_differential_cross_section
from scipy.constants import Avogadro, physical_constants
import time

def plot_spectrum(xx, **params):
    figsize = params.get('figsize', None)
    title = params.get('title', "Proton Energy Spectrum")
    data = xx  # last column
    plt.figure(figsize=figsize)
    plt.hist(data, bins='auto', edgecolor='black')
    # plt.xlabel("proton energy $E_{\\text{p}}\\ [\\text{MeV}]$")
    # plt.ylabel("counts")
    # plt.title(title)
    plt.show()
    plt.close()

class Beam:
    def __init__(self, beam, N_beampart, N_protons, N_neutrons):
        self.list = beam
        self.N_beampart = N_beampart
        self.N_protons = N_protons
        self.N_neutrons = N_neutrons
        self.std_dev_E = None
        self.ave_E = None

    @classmethod
    def generate(cls, target: Target, aperture: Aperture, energy, Npart, cross_section=None,
                 E_scatt=True, E_loss=True):
        '''
        The most realistic way to generate a beam.
        Other methods (useful for testing):
            generate_monoenergetic_parallel,
            generate_monoenergetic_perpendicular.
        '''
        if target == None or aperture == None:
            raise Exception("you CANNOT generate a beam with no Target or Aperture specified!")

        # load cross-section
        if cross_section is None:
            cross_section = ENDF_differential_cross_section.differential_cross_section(
                'E4R84432_e4.endf.endf2gnd.endf')
        f = cross_section.compute_f(energy * 1E6)
        g = cross_section.compute_g(energy * 1E6)
        sigma_s = cross_section.compute_sigma_s(energy * 1E6)

        # generate proton beams at aperture
        protons = np.empty((0, 5))
        N_protons = 0

        # calculate u_min(cos maximum probably scattering angle)
        # u_min = 0.0
        u_min = np.cos(np.arctan((target.max_r + aperture.max_r) / aperture.distance))

        while len(protons) < Npart:
            protons_initial, frac = target.generate_particles(Npart * 5, f, g, u_min)
            start_x, start_y, start_z = protons_initial[:, 0], protons_initial[:, 1], protons_initial[:, 2]
            u, phi = protons_initial[:, 3], protons_initial[:, 4]
            theta = np.arccos(u)
            end_x = start_x + (aperture.distance - start_z) * np.tan(theta) * np.cos(phi)
            end_y = start_y + (aperture.distance - start_z) * np.tan(theta) * np.sin(phi)
            idx = aperture.isPassed(end_x, end_y)
            accepted = protons_initial[idx]

            n_need = Npart - len(protons)
            if len(accepted) >= n_need:
                accepted = accepted[:n_need]
                N_protons += (np.where(idx)[0][n_need - 1] + 1) / frac
                protons = np.vstack([protons, accepted])
                break
            else:
                N_protons += len(protons_initial) / frac
                protons = np.vstack([protons, accepted])

        # calculate neutrons based on the cross-section
        N_neutrons = N_protons * 14.027 / (2.0 * Avogadro) / target.mass_thickness / sigma_s * 1E27

        # calculate proton energy
        start_x, start_y, start_z = protons[:, 0], protons[:, 1], protons[:, 2]
        u, phi = protons[:, 3], protons[:, 4]
        theta = np.arccos(u)
        if E_scatt and E_loss:
            E_par = energy * np.cos(theta) * np.cos(theta)
            E_par -= target.get_ESP(E_par) * ((target.thickness - start_z) * 100) * target.density
        elif E_loss:
            E_par = energy * np.ones_like(theta)
            E_par -= target.get_ESP(E_par) * ((target.thickness - start_z) * 100) * target.density
        elif E_scatt:
            E_par = energy * np.cos(theta) * np.cos(theta)
        else:
            E_par = energy * np.ones_like(theta)

        # generate beam list, [x, a, y, b, t, E]
        beam = np.stack([start_x + (target.thickness - start_z) * np.tan(theta) * np.cos(phi),
                         np.tan(theta) * np.cos(phi),
                         start_y + (target.thickness - start_z) * np.tan(theta) * np.sin(phi),
                         np.tan(theta) * np.sin(phi), np.zeros_like(E_par), E_par], axis=1)  # [x, a, y, b, t, E]
        # std_dev_E = np.std(beam[:, -1])
        # define the class
        obj = cls(beam, len(beam), N_protons, N_neutrons)
        obj.compute_energy_stats()
        # input(f"{obj.N_beampart}, {obj.N_protons}, {obj.N_neutrons}")
        return obj
    @classmethod
    def generate_monoenergetic_parallel(cls, target:Target, energy, Npart):
        '''
        Generate a monoenergetic parallel beam:
        E = E_n, a = b = 0
        '''
        if target == None:
            raise Exception("you CANNOT generate a beam with no Target specified!")
        # generate beam list, [x, a, y, b, t, E]
        beam = np.zeros((Npart, 6))
        beam[:, 0], beam[:, 2] = target.get_initialPosition_N(Npart)
        beam[:, -1] = energy
        obj = cls(beam, len(beam), None, None)
        obj.compute_energy_stats()
        # input(f"{obj.N_beampart}, {obj.N_protons}, {obj.N_neutrons}")
        return obj

    @classmethod
    def generate_monoenergetic_cone(cls, energy, Npart, a_max):
        '''
        Generate a monoenergetic cone-like beam:
        E = E_n, x = y = 0, |a(b)| <= a_max 
        '''
        # generate beam list, [x, a, y, b, t, E]
        beam = np.zeros((Npart, 6))
        beam[:, 1] = np.random.uniform(-a_max, a_max, Npart)
        beam[:, 3] = np.random.uniform(-a_max, a_max, Npart)
        beam[:, -1] = energy
        obj = cls(beam, len(beam), None, None)
        obj.compute_energy_stats()
        # input(f"{obj.N_beampart}, {obj.N_protons}, {obj.N_neutrons}")
        return obj

    def compute_energy_stats(self):
        if self.list is None:
            raise ValueError("Error!! Empty beam list!")
        self.std_dev_E = np.std(self.list[:, -1])
        self.ave_E = np.mean(self.list[:, -1])

    def trans(self, magnet: Magnets):
        """
        Vectorized version of Beam_trans.

        beam: (N, 6) numpy array, columns [x, a, y, b, t, E]
        magnet.TM: (M, 10) numpy array, first 5 columns are coeffs, last 5 columns are powers
        """
        beam = np.asarray(self.list)
        N = beam.shape[0]

        # coordinate = [x, a, y, b, t, delta]
        delta = (beam[:, 5] - magnet.reference_energy) / magnet.reference_energy
        coordinate = np.hstack([beam[:, :5], delta[:, np.newaxis]])  # shape (N,6)

        # TM coeffs and powers
        coeffs = magnet.TM[:, :5]  # shape (M,5)
        powers = magnet.TM[:, 5:]  # shape (M,5)

        # Initialize output
        beam_out = np.zeros((N, 5))
        # Compute monomials
        monomials = np.prod(coordinate[:, None, :] ** powers[None, :, :], axis=2)  # shape (N,M)
        # beam_out: shape (N,5)
        beam_out = monomials @ coeffs

        # switch from (x,a,y,b,t,delta) back to (x,a,y,b,t,E)
        beam_out[:, 5] = (beam_out[:, 5] + 1) * magnet.reference_energy

        beam_new = copy.deepcopy(self)
        beam_new.list = beam_out
        beam_new.compute_energy_stats()
        beam_new.N_beampart = len(beam_out)
        return beam_new
    
    def trans_3d(self, magnet: Magnets, exit_plane: tuple, save_trace = False, save_path = None, batch_size = 1024):
        '''
        Calculate beam trace in 3D Mag-field of magnet using Boris method.
        **Origin must set at the geometric center of the target**
        **(x, y, z) coordinate system, reference beam trace is in xy plane @ z=0**

        magnet: Magnets object with self.type = 3d
        exit_plane: (x0, y0, k) tuple to describe the exit plane after which you should use beam.hit()
        save_trace: if True, save the traces len(beam) particles to save_path.
        batch_size: number of particles to do vectorized calculation.
        '''
        if save_trace == True and save_path is None:
            raise ValueError("save_path of particle's traces is not specified.")
        else:
            if save_trace: trace = []   
            beam = np.asarray(self.list) # coordinate = [x, a, y, b, t, E]
            N = beam.shape[0]
            # transformation from beam reference to (x,y,z,vx,vy,vz)
            beam_in = np.zeros((N, 9)) # [x, y, z, vx, vy, vz, m, t, flag] flag = 1 means arrived
            beam_in[:, 0] = - beam[:, 0]
            beam_in[:, 2] = beam[:, 2]
            m0c2 = physical_constants["proton mass"][0] * physical_constants["speed of light in vacuum"][0] ** 2
            beam_in[:, 6] = (beam[:, -1] * 1e6 * physical_constants["elementary charge"][0] + m0c2) / physical_constants["speed of light in vacuum"][0] ** 2 # mass
            gamma = (beam[:, -1] * 1e6 * physical_constants["elementary charge"][0]+ m0c2) / m0c2
            velocity = np.sqrt(1 - 1/gamma**2) * physical_constants["speed of light in vacuum"][0]
            beam_in[:, 3] = - velocity * beam[:, 1] #vx
            beam_in[:, 5] = velocity * beam[:, 3] #vz
            beam_in[:, 4] = np.sqrt(velocity**2 - beam_in[:, 3]**2 - beam_in[:, 5]**2) #vy
            # calculation
            N_arrived = 0
            beam_out = np.empty((0, 9))
            while N_arrived < N:
                N_batch = min(N - N_arrived, batch_size)
                if save_trace: 
                    beam_batch_out, trace_batch = self._batch_trans_3d(Magnet = magnet, exit_plane = exit_plane, beam_batch_in = beam_in[N_arrived:N_arrived+N_batch], save_trace = save_trace)
                    trace.extend(trace_batch)
                else: 
                    beam_batch_out = self._batch_trans_3d(Magnet = magnet, exit_plane = exit_plane, beam_batch_in = beam_in[N_arrived:N_arrived+N_batch], save_trace = save_trace)
                N_arrived += N_batch
                beam_out = np.vstack((beam_out, beam_batch_out))
            # save_trace
            if save_trace: 
                np.save(os.path.join(save_path, 'traces.npy'), np.array(trace))    
        # transformation from (x,y,z,vx,vy,vz) back to beam reference
        beam_out_ref = np.zeros()
        beam_new = copy.deepcopy(self)
        beam_new.list = beam_out
        beam_new.compute_energy_stats()
        beam_new.N_beampart = len(beam_out)
        return beam_new
    def _batch_trans_3d(self, dt = 2e-11, save_Ninterval = 20, q = physical_constants["elementary charge"][0],**kwargs):
        beam_batch_in = kwargs['beam_batch_in']
        exit_plane = kwargs['exit_plane']
        magnet = kwargs['Magnet']
        save_trace = kwargs['save_trace']
        t_now = 0
        n_iter = 0
        if save_trace: trace = [[np.array(beam_batch_in[i][0], beam_batch_in[i][1], beam_batch_in[i][2], beam_batch_in[i][7])] for i in range(len(beam_batch_in))]
        while np.any(beam_batch_in[:, 8] == 0): # There are still particles not arrived
            # Boris push
            idx = np.where(beam_batch_in[:, 8] == 0)[0] # index of particles not arrived
            B = magnet.get_B_at(np.stack((beam_batch_in[idx, 0], beam_batch_in[idx, 1], beam_batch_in[idx, 2]), axis=1))
            B_magnitude = np.linalg.norm(B, axis=1, keepdims=True)
            mass = beam_batch_in[idx, 6]
            t = np.tan(-q * B_magnitude * dt / mass / 2)
            B_magnitude[B_magnitude == 0] = 1
            b = B / B_magnitude
            M = self._Boris_matrix(t, b, len(idx))
            beam_batch_in[idx, 3:6] = np.einsum('ijk,ik->ij', M, beam_batch_in[idx, 3:6]) # velocity update
            beam_batch_in[idx, 0:3] += beam_batch_in[idx, 3:6] * dt # position update
            beam_batch_in[idx, 7] += dt
            t_now += dt
            n_iter += 1

            # save trace
            if n_iter % save_Ninterval == 0:
                if save_trace:
                    for i in idx:
                        trace[i].append(np.array(beam_batch_in[i][0], beam_batch_in[i][1], beam_batch_in[i][2], beam_batch_in[i][7]))

            # check if particle arrived
            for i in idx:
                if beam_batch_in[i][1] <= exit_plane[2]*(beam_batch_in[i][0] - exit_plane[0]) + exit_plane[1]:
                    beam_batch_in[i][8] = 1
                    if save_trace: trace[i].append(np.array(beam_batch_in[i][0], beam_batch_in[i][1], beam_batch_in[i][2], beam_batch_in[i][7]))
        if save_trace:
            return beam_batch_in, trace
        else:
            return beam_batch_in
    def _Boris_matrix (self, t, b, len):
        M = np.zeros((len, 3, 3))
        C = 2*t/(1+t**2)
        D = 2*t**2/(1+t**2)
        M[:, 0, 0] = 1 + D*(b[:, 0]**2 - 1)
        M[:, 0, 1] = -C*b[:, 2]+D*b[:, 0]*b[:, 1]
        M[:, 0, 2] = C*b[:, 1]+D*b[:, 0]*b[:, 2]
        M[:, 1, 0] = C*b[:, 2]+D*b[:, 0]*b[:, 1]
        M[:, 1, 1] = 1 + D*(b[:, 1]**2 - 1)
        M[:, 1, 2] = -C*b[:, 0]+D*b[:, 1]*b[:, 2]
        M[:, 2, 0] = -C*b[:, 1]+D*b[:, 0]*b[:, 2]
        M[:, 2, 1] = C*b[:, 0]+D*b[:, 1]*b[:, 2]
        M[:, 2, 2] = 1 + D*(b[:, 2]**2 - 1)
        return M

    def plot_energy_spectrum(self, **params):
        bin_width = params.get('bin_width')
        figsize = params.get('figsize', None)
        title = params.get('title', "Proton Energy Spectrum")
        filename = params.get('filename')
        data_filename = params.get('data_filename')
        show = params.get('show', False)
        if self.list is None:
            raise ValueError("Beam.list is empty")

        data = self.list[:, -1]
        plt.figure(figsize=figsize)
        if bin_width:
            bin_start = np.floor(data.min() / bin_width) * bin_width
            bin_end = np.ceil(data.max() / bin_width) * bin_width
            bins = np.arange(bin_start, bin_end + bin_width, bin_width)
        else:
            bins = 'auto'
        counts, bin_edges, _ = plt.hist(data, bins=bins, edgecolor='black')
        plt.xlabel("proton energy $E_{\\text{p}}\\ [\\text{MeV}]$")
        plt.ylabel("counts")
        plt.title(title)
        plt.tight_layout()
        if filename:
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            plt.savefig(f"{filename}.png", dpi=300, bbox_inches='tight')
        if show:
            plt.show()
        plt.close()
        if data_filename:
            os.makedirs(os.path.dirname(data_filename), exist_ok=True)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
            data_output = np.column_stack((bin_centers, counts))
            np.savetxt(data_filename, data_output, delimiter='\t', fmt=("%.3f", "%.1f"))

    def hit(self, focalplane: Focalplane):
        # shape: [N_part, N_node], directed length of nodes to particle trajectories
        node_distances_x = (focalplane.geometry[None, :, 0] - self.list[:, None, 0] -
                            self.list[:, None, 1] * (focalplane.geometry[None, :, 1]))

        # determine the number of intersection points
        mask = node_distances_x >= 0
        transitions = np.diff(mask.astype(int), axis=1)  # [N_part, N_node-1], +1/-1 means the trajectory pass through the two nodes
        has_transition = transitions != 0
        # transition_indices = [np.where(row != 0)[0] for row in transitions]
        N_intersection_points = np.sum(has_transition, axis=1)
        if np.any(N_intersection_points >= 2):
            raise ValueError("Too many transitions (>=2)!")
        # For cases where only one intersection is allowed, find the index position
        # Method: argmax returns the first True position, but for all False it gives 0, which needs to be masked (set to -1)
        first_transition_idx = np.argmax(has_transition, axis=1)
        first_transition_idx_selected = first_transition_idx[np.any(has_transition, axis=1)]
        particle_selected = self.list[np.any(has_transition, axis=1)]
        node_distance_x_selected = node_distances_x[np.any(has_transition, axis=1), :]  # [N_hits, N_node]
        N_hits = np.count_nonzero(np.any(has_transition, axis=1))
        rows = np.arange(N_hits)

        d0 = node_distance_x_selected[rows, first_transition_idx_selected]
        d1 = node_distance_x_selected[rows, first_transition_idx_selected + 1]

        g0 = focalplane.geometry[first_transition_idx_selected]
        g1 = focalplane.geometry[first_transition_idx_selected + 1]
 
        x_hits, z_hits, l_hits = ((-d0[:, None] * g1 + d1[:, None] * g0) /
                                  (-d0 + d1)[:, None]).T
        y_hits = particle_selected[:, 2] + particle_selected[:, 3] * z_hits

        return ProtonHits(x_hits, y_hits, z_hits, l_hits)

class ProtonHits:
    def __init__(self, x_hits, y_hits, z_hits, l_hits):
        self.x_hits, self.y_hits, self.z_hits = x_hits, y_hits, z_hits
        self.l_hits = l_hits
        self.N_hits = len(self.l_hits)
        self.std_dev_l = np.std(self.l_hits)
        self.x_mean, self.y_mean, self.z_mean = np.mean(self.x_hits), np.mean(self.y_hits), np.mean(self.z_hits)
        self.l_mean = np.mean(self.l_hits)
