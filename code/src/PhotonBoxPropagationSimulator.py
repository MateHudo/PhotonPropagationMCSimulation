import numpy as np
from datetime import datetime as dt
import matplotlib.pyplot as plt

class photon_box_propagation_simulator:
        """
    Simulates photon propagation in a 3D box-shaped medium, tracking each step,
    considers photoeffect (PE), Compton scattering (CS), and pair production (PP). Assumes local
    electron energy deposition and no braking radiation or fluorescence photons (for the latter
    fluorescence yield is important. In low-Z matter like tissue it is ok to neglect it, but in
    high-Z material like Lead, fluo. yield is high, but those photons are be absorbed locally). 
    Data for attenuation coefficients is gotten from NIST XCOM, using soft tissue mixture.
    todo: incorporate NIST raw data and use LACloader method to calculate LAC from MAC!
    """
    def __init__(self, BoxData, LAC_data, Emin_terminate=0.001):
        self.BoxData = BoxData
        self.LAC_data = LAC_data
        self.Emin_terminate = Emin_terminate
        
        # Define constants once at initialization instead of in methods
        self.m_e = 0.511  # Electron rest mass energy in MeV
        self.pair_production_threshold = 1.022  # MeV
        self.annihilation_energy = 0.511  # MeV per gamma in pair production
        
        # Pre-extract frequently used LAC data for efficiency
        self.energy_list = self.LAC_data['energy']
        self.lac_cs_list = self.LAC_data['lac_cs']
        self.lac_pe_list = self.LAC_data['lac_pe'] 
        self.lac_pp_list = self.LAC_data['lac_pp']
        
        # Pre-extract box geometry data
        self.cuboid_dimensions = self.BoxData['cuboid_dimensions']
        self.face_normals = self.BoxData['face_normals']
        self.r_faces = self.BoxData['r_faces']
        self.n_hvl = self.BoxData['n_hvl']  # Number of HVLs for the simulation
        self.n_hvl_leakage = self.BoxData['n_hvl_leakage']

    def SphericalToCartesian(self, theta, phi):
        x = np.cos(theta)
        y = np.sin(theta) * np.cos(phi)
        z = np.sin(theta) * np.sin(phi)
        return np.array([x, y, z])

    def NewDirectionUnitVector(self, u_old, theta_sc, phi_sc):
        u_old = u_old / np.linalg.norm(u_old) # necessary?
        u_sc_local = self.SphericalToCartesian(theta_sc, phi_sc)
        # if direction before scattering is close to initial direction of the beam, simple transform
        if np.allclose(u_old, [1, 0, 0]):
            return u_sc_local
        ref_vec = np.array([1, 0, 0])
        v = np.cross(ref_vec, u_old)
        v /= np.linalg.norm(v) # necessary?
        w = np.cross(u_old, v)
        R = np.column_stack((u_old, v, w))
        u_new = R @ u_sc_local
        u_new /= np.linalg.norm(u_new)
        return u_new

    def ComptonScatteringInteraction(self, E, u_particle):
        # Use pre-defined constant instead of defining it each time
        alpha = E / self.m_e
        eps_min = 1 / (1 + 2 * alpha)
        beta1 = - np.log(eps_min)
        beta2 = 0.5 * (1 - eps_min ** 2)
        while True:
            r1, r2, r3 = np.random.uniform(0, 1, 3)
            if r1 < beta1 / (beta1 + beta2):
                eps_prop = eps_min ** r2
            else:
                eps_prop = np.sqrt(eps_min ** 2 + (1 - eps_min ** 2) * r2)
            t = (1 - eps_prop) / alpha / eps_prop
            g_eps = 1 - eps_prop * t * (2 - t) / (1 + eps_prop ** 2)
            if g_eps > r3:
                eps = eps_prop
                break
        E_scattered = E * eps
        theta_cs = np.arccos(1 - 1 / alpha * (1 / eps - 1))
        phi_cs = 2 * np.pi * np.random.uniform(0, 1)
        u_new = self.NewDirectionUnitVector(u_particle, theta_cs, phi_cs)
        return E_scattered, u_new

    def PhotonFreePath(self, E):
        # Use pre-extracted data instead of accessing dictionary each time
        lac_cs = np.interp(E, self.energy_list, self.lac_cs_list)
        lac_pe = np.interp(E, self.energy_list, self.lac_pe_list)
        u1, u2, u3 = np.random.uniform(0, 1, 3)
        s_cs = - np.log(u1) / lac_cs
        s_pe = - np.log(u2) / lac_pe
        if E > self.pair_production_threshold:  # Use pre-defined constant
            lac_pp = np.interp(E, self.energy_list, self.lac_pp_list)
            s_pp = - np.log(u3) / lac_pp
        else:
            s_pp = np.inf
        free_path_length = min(s_cs, s_pe, s_pp)
        if free_path_length == s_cs:
            interaction_type = 'compton'
        elif free_path_length == s_pe:
            interaction_type = 'phot'
        else:
            interaction_type = 'pair'
        return free_path_length, interaction_type

    def BoxPathLength3D(self, r_particle, u_particle):
        min_path_length = np.inf
        exit_plane_index = -1
        for plane_index in range(6):
            r_plane, n_plane = self.r_faces[plane_index], self.face_normals[plane_index]
            relative_direction_scaling_factor = np.dot(u_particle, n_plane)
            if relative_direction_scaling_factor <= 0:
                continue
            path_length_to_plane = np.abs(np.dot(n_plane, r_particle - r_plane)) / relative_direction_scaling_factor
            if path_length_to_plane < min_path_length:
                min_path_length = path_length_to_plane
                exit_plane_index = plane_index
        return min_path_length, exit_plane_index

    def ParticleStep(self, PreStepParticleInfo):
        r_particle = PreStepParticleInfo['r']
        u_particle = PreStepParticleInfo['u']
        E = PreStepParticleInfo['E']
        eid = PreStepParticleInfo['eid']
        tid = PreStepParticleInfo['tid']

        add_photons_to_simulate = [] 
        box_path_length, exit_plane_index = self.BoxPathLength3D(r_particle, u_particle)
        free_path_length, interaction_type = self.PhotonFreePath(E)
        
        if free_path_length > box_path_length:
            interaction_type = 'exit'
            step_length = box_path_length
            E_dep = 0.0
            E_new = E
            #print(f"   --->Particle exited the box at face index {exit_plane_index} and energy {E_new:.4f} MeV")
            u_new = u_particle
            r_new = r_particle + box_path_length * u_particle
            particle_termination = True
            #event_termination = True  ## this is more complex, we will set it later
        else:
            r_new = r_particle + free_path_length * u_particle
            step_length = free_path_length
            exit_plane_index = -1
            if interaction_type == 'compton':
                E_new, u_new = self.ComptonScatteringInteraction(E, u_particle)
                E_dep = E - E_new
                # only in CS, the condition about termination based on energy is checked, is E is low enough, we terminate the particle
                if E_new < self.Emin_terminate:
                    particle_termination = True
                # if E_new >= self.Emin_terminate; then we add the photon to the list to simulate
                else:
                    particle_termination = False
                    add_photons_to_simulate.append({
                        'eid': eid,
                        'tid': tid,
                        'E': E_new,
                        'u': u_new,
                        'r': r_new,
                    })
            elif interaction_type == 'phot':
                E_new = 0.0
                E_dep = E
                u_new = None
                particle_termination = True
            elif interaction_type == 'pair':
                E_new = 0.0
                E_dep = E - self.pair_production_threshold  # Use pre-defined constant
                u_new = None
                particle_termination = True
                phi = np.random.uniform(0, 2 * np.pi)
                th = np.arccos(np.random.uniform(-1, 1))
                x, y, z = np.sin(th) * np.cos(phi), np.sin(th) * np.sin(phi), np.cos(th)
                u_pp = np.array([x, y, z])
                add_photons_to_simulate.extend([
                    {'eid': eid, 'tid': tid + 1, 'E': self.annihilation_energy, 'u': u_pp, 'r': r_new},  # Use pre-defined constant
                    {'eid': eid, 'tid': tid + 2, 'E': self.annihilation_energy, 'u': -u_pp, 'r': r_new}  # Use pre-defined constant
                ])
            elif interaction_type == 'exit':
                pass  # already set above
            else:
                raise RuntimeError("Unknown interaction type")
            
        PostStepParticleInfo = {
            'eid': eid,
            'tid': tid,
            'r': r_new,
            'u': u_new,
            'E': E_new,
        }
        
        return {
            'eid': eid,
            'tid': tid,
            
            #'PreStepParticleInfo': PreStepParticleInfo,  # Original particle info before the step
            #'PostStepParticleInfo': PostStepParticleInfo,  # Particle info after the step

            'E_prestep': E,  # Original energy before the step
            'r_prestep': r_particle,  # Original position before the step
            'u_prestep': u_particle,  # Original direction before the step
            
            #'free_path_length': free_path_length,
            #'box_path_length': box_path_length,
            #'step_length': step_length,
            
            'interaction': interaction_type,
            'exit_plane_index': exit_plane_index,
            'E_dep': E_dep,
            'add_photons_to_simulate': add_photons_to_simulate,

            'E_new': E_new,
            'u_new': u_new,
            'r_new': r_new,
            

            #'particle_termination': particle_termination,
            #'event_termination': event_termination,  ## this is not correct !!!!
            
            
            # anything else???
        }

    def PhotonPropagationSimulation(self, Nsim, SourcePhotons_template):
        """ Simulates the propagation of photons through a defined box geometry.
        Args:
            Nsim (int): Number of simulation events.
            SourcePhotons_template (dict): Template for the source photons, containing 'r', 'u', and 'E'.
        """
        Ephotons = SourcePhotons_template['E']
        #n_hvl = self.n_hvl
        print(f"Starting simulation: E={Ephotons}; n_hvl = {self.n_hvl}; {Nsim} events\n\n")

        # create a timestamped filename for saving simulation info
        timestamp_str = dt.now().strftime("%Y%m%d_%H%M%S")
        addinfo_filename = f"output_data/run_{timestamp_str}.txt"
        
        # lets add track of total energy that exists on the right size of the box - build-up contribution
        self.E_out_tot = 0.0
        self.E_deposited_tot = 0.0
        self.E_backscattered = 0.0
        self.E_leakage = 0.0

        for eid in range(Nsim):
            # Print progress every 10% of Nsim
            if eid % (Nsim // 10) == 0:
                print(f"Progress: {eid / Nsim * 100:.2f}%") 

            photons_to_simulate = [{
                'eid': eid,
                'tid': 0,
                'r': SourcePhotons_template['r'],
                'u': SourcePhotons_template['u'],
                'E': SourcePhotons_template['E'],
            }]
            
            
            while photons_to_simulate:
                # Process the first photon in the list
                PreStepData = photons_to_simulate.pop(0)
                # simulate the particle step
                result = self.ParticleStep(PreStepData)
                # extend the list of photons to simulate with new photons generated in this step
                photons_to_simulate.extend(result['add_photons_to_simulate'])
                
                # extract some data from the photon step 
                exit_plane_index = result['exit_plane_index']
                Edep = result['E_dep']
                Enew = result['E_new']



                #print(f"Event {result['eid']} | Track {result['tid']} | Interaction: {result['interaction']} | E_dep: {result['E_dep']:.4f} MeV")
                #print(result)
                #print("-->",len(photons_to_simulate), "photons to simulate after this step")
                

                # write the result to the file
                with open(addinfo_filename, "a") as f:
                    
                    # i) A more detailed output
                    # xint,yint,zint = result['r_new']
                    # num_photons_to_simulate = len(photons_to_simulate)
                    # f.write(
                    #     f" {result['eid']} {result['tid']} {result['interaction']} {exit_plane_index} "
                    #     f"{result['E_dep']:.4f} {result['E_new']} {xint} {yint} {zint} "
                    #     f"{result['u_new']} {num_photons_to_simulate}\n"
                    # )
                
                    # ii) A more compact output
                    xint, yint, zint = result['r_new']
                    # # this needs change
                    # f.write(
                    #     f" {result['eid']} {result['tid']} {result['interaction']} {result['E_dep']:.4f}"
                    #     f" {xint} {yint} {zint}\n"    
                    # )

                    # iii) Output for Edep distribution
                    xint, yint, zint = result['r_new']
                    f.write(
                        f"{Edep} {xint} {yint} {zint}\n"
                    )

                    ## iv) ...


                ## Build-up contribution, also add other faces and deposited energy
                # a) energy deposited insude the box
                if exit_plane_index == -1:
                    self.E_deposited_tot += Edep

                # b) - build-up contribution
                elif exit_plane_index == 0:  # Assuming exit_plane_index 0 is the right side of the box
                    self.E_out_tot += Enew
                    #print(f"Total energy on the right side of the box: {self.E_out_tot:.4f} MeV")

                # c) backscattered photons - exit on the face where they entered
                elif exit_plane_index == 1:  # Assuming exit_plane_index 2 is the back side of the box
                    self.E_backscattered += Enew
                    
                # d) leakage photons - exit on the side faces of the box
                else:
                    self.E_leakage += Enew

