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
    def __init__(self, box_num_hvls, lac_loader, GPS_template, Emin_terminate=0.001):
        # number of HVLs in x, y, z directions
        self.num_hvlX, self.num_hvlY, self.num_hvlZ = box_num_hvls
        # linear attenuation coefficients loader class instance
        self.lac_loader = lac_loader
        self.lac_energy, self.lac_pp, self.lac_cs, self.lac_pe, self.lac_total = lac_loader.get_lac_data()
        self.Emin_terminate = Emin_terminate
        
        # Define constants once at initialization instead of in methods
        self.m_e = 0.511  # Electron rest mass energy in MeV
        self.pair_production_threshold = 1.022  # MeV
        self.annihilation_energy = 0.511  # MeV per gamma in pair production

        # GPS_template should contain:
        # "E" (MeV); "r" -  (), "u" ....
        # GPS = General Particle (Photon) Source - monochromatic pencil beam for now ...
        self.GPS_template = GPS_template
        self.E0 = GPS_template['E0']
        self.r_entrance = GPS_template['r_entrance']
        self.u_entrance = GPS_template['u_entrance']

        # Define box geometry and related parameters based on inputs pars: E0 and n_hvls
        self.define_box_related_parameters()


        # manual calculation of latter 2 parametrs!
        self.r_entrance = np.array([-self.box_sizeX/2 + 1e-10,0,0])
        self.u_entrance = np.array([1,0,0])
        print(f"GPS_template: \nE0={self.E0}; r_entr={self.r_entrance}; u_entr={self.u_entrance}")

    def define_box_related_parameters(self):
        self.gps_hvl = np.log(2) / np.interp(self.E0, self.lac_energy, self.lac_total)  # in cm
        self.box_sizeX = self.num_hvlX * self.gps_hvl  # cm
        self.box_sizeY = self.num_hvlY * self.gps_hvl  # cm
        self.box_sizeZ = self.num_hvlZ * self.gps_hvl  # cm
        #self.BoxData = self.BoxGeometry()

        #* face normals 
        nF = np.array([1, 0, 0])  # normal vector of front face
        nB = np.array([-1, 0, 0])  # normal vector of back face
        nL = np.array([0, 1, 0])  # normal vector of left face
        nR = np.array([0, -1, 0])  # normal vector of right face
        nT = np.array([0, 0, 1])  # normal vector of top face
        nBo = np.array([0, 0, -1])  # normal vector of bottom face
        self.face_normals = np.array([nF, nB, nL, nR, nT, nBo])

        #* faces positions - center points on the planes
        rF = np.array([self.box_sizeX / 2, 0, 0])  # Front face position
        rB = np.array([-self.box_sizeX / 2, 0, 0]) # Back face position
        rL = np.array([0, self.box_sizeY / 2, 0])  # Left face position
        rR = np.array([0, -self.box_sizeY / 2, 0])  # Right face position
        rT = np.array([0, 0, self.box_sizeZ / 2])  # Top face position
        rBo = np.array([0, 0, -self.box_sizeZ / 2])  # Bottom face position
        self.face_centers = np.array([rF, rB, rL, rR, rT, rBo])

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
        lac_cs = np.interp(E, self.lac_energy, self.lac_cs)
        lac_pe = np.interp(E, self.lac_energy, self.lac_pe)
        u1, u2, u3 = np.random.uniform(0, 1, 3)
        s_cs = - np.log(u1) / lac_cs
        s_pe = - np.log(u2) / lac_pe
        if E > self.pair_production_threshold:  # Use pre-defined constant
            lac_pp = np.interp(E, self.lac_energy, self.lac_pp)
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
            r_plane, n_plane = self.face_centers[plane_index], self.face_normals[plane_index]
            #print(f"plane_index:{plane_index}; u_particle:{u_particle}; n_plane={n_plane}")
            relative_direction_scaling_factor = np.dot(u_particle, n_plane)
            #print(f"re_dir_sc_factor:{relative_direction_scaling_factor}")
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
        # fixme: could remove eid and tid for simple simulation of build-up!!
        eid = PreStepParticleInfo['eid']
        tid = PreStepParticleInfo['tid']

        add_photons_to_simulate = [] 
        box_path_length, exit_plane_index = self.BoxPathLength3D(r_particle, u_particle)
        free_path_length, interaction_type = self.PhotonFreePath(E)
        
        if free_path_length > box_path_length:
            interaction_type = 'exit'
            #step_length = box_path_length
            E_dep = 0.0
            E_new = E
            #print(f"   --->Particle exited the box at face index {exit_plane_index} and energy {E_new:.4f} MeV")
            u_new = u_particle
            r_new = r_particle + box_path_length * u_particle
            #particle_termination = True
            #event_termination = True  ## this is more complex, we will set it later
        else:
            r_new = r_particle + free_path_length * u_particle
            #step_length = free_path_length
            exit_plane_index = -1
            if interaction_type == 'compton':
                E_new, u_new = self.ComptonScatteringInteraction(E, u_particle)
                E_dep = E - E_new
                # only in CS, the condition about termination based on energy is checked, is E is low enough, we terminate the particle
                if E_new < self.Emin_terminate:
                    particle_termination = True
                # if E_new >= self.Emin_terminate; then we add the photon to the list to simulate
                else:
                    #particle_termination = False
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
                #particle_termination = True
            elif interaction_type == 'pair':
                E_new = 0.0
                E_dep = E - self.pair_production_threshold  # Use pre-defined constant
                u_new = None
                #particle_termination = True
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
            
        # PostStepParticleInfo = {
        #     'eid': eid,
        #     'tid': tid,
        #     'r': r_new,
        #     'u': u_new,
        #     'E': E_new,
        # }
        
        return {
            #* main info
            'interaction': interaction_type,
            'exit_plane_index': exit_plane_index,
            'E_dep': E_dep,
            'add_photons_to_simulate': add_photons_to_simulate,
            'E_new': E_new,
            'u_new': u_new,
            'r_new': r_new,
            
            #* "extras - maybe useful one day"
            # 'eid': eid,
            # 'tid': tid,
            # 'E_prestep': E,  # Original energy before the step
            # 'r_prestep': r_particle,  # Original position before the step
            # 'u_prestep': u_particle,  # Original direction before the step           
 
            #'PreStepParticleInfo': PreStepParticleInfo,  # Original particle info before the step
            #'PostStepParticleInfo': PostStepParticleInfo,  # Particle info after the step
            #'free_path_length': free_path_length,
            #'box_path_length': box_path_length,
            #'step_length': step_length,
            #'particle_termination': particle_termination,
            #'event_termination': event_termination,  ## this is not correct !!!!
            #'interaction_type': interaction_type,
            #'exit_plane_index': exit_plane_index,
            # anything else???
        }

    def simulate(self, Nsim,doWrite=False):
        """ Simulates the propagation of photons through a defined box geometry.
        Args:
            Nsim (int): Number of simulation events.
            SourcePhotons_template (dict): Template for the source photons, containing:
                - 'E': Energy of the photons (MeV).
                - 'r': Initial position of the photons (for now center left plane: [-a/2,0,0] ~ pencil-beam)
                - 'u': Initial direction unit vector of the photons (for now along the x-axis: [1,0,0]).
        """
        self.Nsim = Nsim
        print(f"Starting simulation: E={self.E0}; num_hvlX = {self.num_hvlX}; {self.Nsim} events\n\n")
        
        # create a timestamped filename for saving simulation info
        timestamp_str = dt.now().strftime("%Y%m%d_%H%M%S")
        output_file= f"initial_output/run_{timestamp_str}.txt"
        #define a writer function depending on doWrite
        if doWrite:
            def write_line(msg):
                with open(output_file,"a") as f:
                    f.write(msg)
        else:
            def write_line(msg):
                pass # no-op

        # lets add track of total energy that exists on the right size of the box - build-up contribution
        self.E_out_tot = 0.0
        # energy, deposited inside the box
        self.E_deposited_tot = 0.0
        # backscattered energy - exit on the face where they entered
        self.E_backscattered = 0.0
        # leakage energy - exit on the side faces of the box - unwanted i guess
        self.E_leakage = 0.0


        # define printing progress frequency
        progress_report_freq = 10  # %
        eid_progress_set = set(
            int(Nsim * pct / 100)
            for pct in range(progress_report_freq, 101, progress_report_freq)
        )

        #* loop over all simulated photons
        for eid in range(Nsim):
            # Print progress every {progress_report_freq}% of Nsim
            if eid in eid_progress_set:
                print(f"\tProgress: {eid / Nsim * 100:.2f}%")

            # add initial photon to simulate
            photons_to_simulate = [{
                'eid': eid,
                'tid': 0,
                'r': self.r_entrance,
                'u': self.u_entrance,
                'E': self.E0,
                # fixme: problem with calculating r_entrance
                # 'r': self.GPS_template['r_entrance'],
                # 'u': self.GPS_template['u_entrance'],
                # 'E': self.GPS_template['E0'],
            }]
            

            #* loop that continues until photons_to_simulate is empty
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

                #* optional print of the step result
                #print(f"Event {result['eid']} | Track {result['tid']} | Interaction: {result['interaction']} | E_dep: {result['E_dep']:.4f} MeV")
                #print(result)
                #print("-->",len(photons_to_simulate), "photons to simulate after this step")
                

                # write the result to the file
                # i) short output
                xint, yint, zint = result['r_new']
                int_type = result['interaction']
                write_line(f"{int_type} {Edep} {xint} {yint} {zint}\n")
                # ii) detailed output
                # num_photons_to_simulate = len(photons_to_simulate)
                # write_line(
                #     f" {result['eid']} {result['tid']} {result['interaction']} {exit_plane_index} "
                #     f"{result['E_dep']:.4f} {result['E_new']} {xint} {yint} {zint} "
                #     #f"{result['u_new']} {num_photons_to_simulate}\n"
                # )
                #iii) ...


                ##* Build-up contribution, also add other faces and deposited energy
                # a) energy deposited inside the box
                if exit_plane_index == -1:
                    self.E_deposited_tot += Edep
                # b) - build-up contribution
                elif exit_plane_index == 0:  # Assuming exit_plane_index 0 is the right side of the box
                    self.E_out_tot += Enew

                # c) backscattered photons - exit on the face where they entered
                elif exit_plane_index == 1:  # Assuming exit_plane_index 2 is the back side of the box
                    self.E_backscattered += Enew
                    
                # d) leakage photons - exit on the side faces of the box
                else:
                    self.E_leakage += Enew
        
        print("Simulation completed.\n")
        # final report of the results
        self.report_results()

    def report_results(self):
        """
        Print out simultion results.
        """
        if not hasattr(self, 'E_deposited_tot'):
            raise RuntimeError("Simulation not yet run. Please run simulate() before reporting results.")

        print("Simulation Results:")
        print(f"  Deposited Energy:   {self.E_deposited_tot:.6f} MeV")
        print(f"  Backscattered:      {self.E_backscattered:.6f} MeV")
        print(f"  Out Total:          {self.E_out_tot:.6f} MeV")
        print(f"  Leakage:            {self.E_leakage:.6f} MeV")
        # add checking whether all energy was counted in simulaion process
        E_tot_simulated = self.Nsim * self.E0
        E_tot_detected = self.E_deposited_tot + self.E_backscattered + self.E_out_tot + self.E_leakage
        print(f"  Total Simulated:    {E_tot_simulated:.6f} MeV")
        print(f"  Total Detected:     {E_tot_detected:.6f} MeV")
        print(f"  Energy Balance:     {E_tot_simulated - E_tot_detected:.6f} MeV (should be close to 0)")
        # calculate buildup
        Eout_theory = E_tot_simulated * 2**(-self.num_hvlX)
        buildup = self.E_out_tot / Eout_theory
        print(f"  Theoretical Out:    {Eout_theory:.6f} MeV")

        # report buildup factor with ifno about E0, n_hvlX, Nsim
        print(f"\nBuildup Factor Report:")
        print(f"  Incident Energy:    {self.E0:.4f} MeV")
        print(f"  Number of HVLs:     {self.num_hvlX}")
        print(f"  Number of Events:   {self.Nsim}")
        print(f"  Buildup Factor:     {buildup:.6f}")