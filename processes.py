import numpy as np
from components import Target, Aperture, Magnets, Focalplane

def Beam_init(method: str, target = None, aperture = None, E_loss = True):
    if method == 'import':
        pass

    if method == 'generate':
        if target == None or aperture == None:
            raise Exception("you CANNOT generate a beam with no Target or Aperture specified!")
        else:
            print("input minimum n-beam energy: [MeV]")
            Emin = float(input())
            print("input maximum n-beam energy: [MeV]")
            Emax = float(input())
            print("input number of bins: ")
            Nbins = int(input())
            print("input number of particles per bin: ")
            Npart = int(input())

            # starting generating the beam
            Beam = []
            E_lst = np.linspace(Emin, Emax, Nbins) # [MeV]
            for E in E_lst:
                N_bin_valid = 0
                while (N_bin_valid < Npart):
                    # random starting position in the target
                    start_r = np.random.uniform(0, target.r)
                    start_phi = np.random.uniform(0, 2*np.pi)
                    start_x = start_r * np.cos(start_phi)
                    start_y = start_r * np.sin(start_phi)
                    start_z = np.random.uniform(0, target.thickness)

                    # random scattering angle
                    u = np.random.uniform(0, 1)
                    theta = np.arcsin(np.sqrt(u))
                    phi = np.random.uniform(0, 2*np.pi)

                    # calculate final position
                    end_x = start_x + (aperture.distance - start_z)*np.tan(theta)*np.cos(phi)
                    end_y = start_y + (aperture.distance - start_z)*np.tan(theta)*np.sin(phi)

                    # check if the particle is inside the aperture
                    if (end_x**2 + end_y**2 > aperture.r**2):
                        print(f"r = {np.sqrt(end_x**2 + end_y**2):.2f} m, particle is outside the aperture")
                        continue
                    else:
                        # calculate energy loss
                        if E_loss:
                            E_par = E * np.cos(theta)
                            E_par -= target.get_ESP(E_par) * ((target.thickness - start_z) * 100) * target.density
                        else:
                            E_par = E

                        # add the particle to the beam
                        Beam.append([start_x + (target.thickness - start_z)*np.tan(theta)*np.cos(phi), np.tan(theta)*np.cos(phi), start_y + (target.thickness - start_z)*np.tan(theta)*np.sin(phi), np.tan(theta)*np.sin(phi), 0, E_par]) #[x, a, y, b, t, E]
                        N_bin_valid += 1
                        print(f"One particle passed! N_bin_valid = {N_bin_valid}, E = {E_par:.2f} MeV")

            return Beam
        
def Beam_trans (magnet: Magnets, beam: list):
    beam_out = []
    for particle in beam:
        coordinate = particle[:5]
        coordinate.append((particle[-1] -  magnet.reference_energy)/magnet.reference_energy) # [x, a, y, b, t, delta]

        particle_out = [0, 0, 0, 0, 0] # [x, a, y, b, t]
        for line in magnet.TM:
            for i, coeff in enumerate(line[:5]):
                if coeff != 0:
                    temp = 1
                    for (coor, power) in zip(coordinate, line[5:]):
                        temp *= coor**power
                    particle_out[i] += temp * coeff
        particle_out.append(particle[-1]) # [x, a, y, b, t, E] at Exit
        beam_out.append(particle_out)
    return beam_out

def Aline (x, y, A, B, C):
    return A*x + B*y + C
def Beam_hit (magnet: Magnets, focalplane: Focalplane, beam: list):
    record = [] # list of (l, y) in [m]
    count = 0

    for particle in beam:
        for (node1, node2) in zip(focalplane.geometry[:-1], focalplane.geometry[1:]):
            if Aline(node1[0], node1[1], -1, particle[1], particle[0] - particle[1]*magnet.length) == 0:
                record.append([node1[2], particle[2] + (node1[1] - magnet.length) * particle[3]])
                count += 1
                break
            elif Aline(node1[0], node1[1], -1, particle[1], particle[0] - particle[1]*magnet.length) * Aline(node2[0], node2[1], -1, particle[1], particle[0] - particle[1]*magnet.length) < 0:
                # calculate intersection point
                denominator = particle[1]*(node2[1] - node1[1])/(node2[0] - node1[0]) - 1
                x0 = (particle[1] * (node1[0] * (node2[1] - node1[1]) / (node2[0] - node1[0]) - node1[1]) - (particle[0] - particle[1]*magnet.length)) / denominator
                z0 = ((node1[0] * (node2[1] - node1[1]) / (node2[0] - node1[0]) - node1[1]) - (node2[1] - node1[1]) / (node2[0] - node1[0]) * (particle[0] - particle[1]*magnet.length)) / denominator

                record.append([node1[2] + np.sqrt((x0 - node1[0]) ** 2 + (z0 - node1[1])**2), particle[2] + particle[3] * (z0 - magnet.length)])
                count += 1
                break

    print(f"{count}/{len(beam)} particles hit the focal plane.")

    return record