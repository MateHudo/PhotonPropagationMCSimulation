import numpy as np
from datetime import datetime as dt
import matplotlib.pyplot as plt

class photon_box_propagation_simulator_buildup_set:
    """
    Simulates photon propagation in a 3D box-shaped medium, tracking each step,
    considers photoeffect (PE), Compton scattering (CS), and pair production (PP). Assumes local
    electron energy deposition and no braking radiation or fluorescence photons (for the latter
    fluorescence yield is important. In low-Z matter like tissue it is ok to neglect it, but in
    high-Z material like Lead, fluorescence yield is high, but those photons are be absorbed locally). 
    Data for attenuation coefficients is gotten from: 
    NIST XCOM (https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html);
    using soft tissue mixture.
    * AIM: to calculate the build-up factor for photons propagating through a box-shaped medium, 
    * and to understand the contribution of different interaction types (PE, CS, PP) to the 
    * energy deposition and photon propagation.
    #! This code optimizes for buildup calculation. For more general simulation code look in PhotonBoxPropagationSimulator.py
    """
    def __init__(self, box_num_hvls, lac_loader, E0, Emin_terminate=0.001):
        # number of HVLs in x, y, z directions
        self.num_hvlX, self.num_hvlY, self.num_hvlZ = box_num_hvls
        # linear attenuation coefficients loader class instance
        self.lac_loader = lac_loader
        self.lac_energy, self.lac_pp, self.lac_cs, self.lac_pe, self.lac_total = lac_loader.get_lac_data()
        self.Emin_terminate = Emin_terminate # minimal energy of photons to terminate simulation (unit: MeV)
        self.E0 = E0

        # Define constants once at initialization instead of in methods
        self.m_e = 0.511  # Electron rest mass energy in MeV
        self.pair_production_threshold = 1.022  # MeV
        self.annihilation_energy = 0.511  # MeV per gamma in pair production

        # Define box geometry and related parameters based on inputs pars: E0 and n_hvls
        self.define_box_related_parameters()
        # manual calculation of latter 2 parameters!
        self.r_entrance = np.array([-self.box_sizeX/2 + 1e-10,0,0])
        self.u_entrance = np.array([1,0,0])
        
        # GPS = General Particle (Photon) Source - monochromatic pencil beam entering at left box surface for now ...
        # GPS_template should contain:
        # * "E0" - initial photon energy (MeV);
        # * "r_entrance" -  entrance position (cm), 
        # * "u_entrance" - unit vector of photon direction
        # ...
        self.GPS_template = {
            'E0': self.E0,
            'r_entrance': self.r_entrance,
            'u_entrance': self.u_entrance,
        }   
        print(f"GPS_template: \n\tE0={self.E0}; \n\tr_entr={self.r_entrance}; \n\tu_entr={self.u_entrance}")
        print("Box HVL (half-value layer) = {:.4f} cm".format(self.gps_hvl))

    def define_box_related_parameters(self):
        self.gps_hvl = np.log(2) / np.interp(self.E0, self.lac_energy, self.lac_total)  # in cm
        self.box_sizeX = self.num_hvlX * self.gps_hvl  # cm
        self.box_sizeY = self.num_hvlY * self.gps_hvl  # cm
        self.box_sizeZ = self.num_hvlZ * self.gps_hvl  # cm

        #* face normals - with indexing order: front, back, left, right, top, bottom [IMPORTANT]
        self.nF = np.array([1, 0, 0])  # normal vector of front face  - indexed as 0
        self.nB = np.array([-1, 0, 0])  # normal vector of back face - indexed as 1
        self.nL = np.array([0, 1, 0])  # normal vector of left face - indexed as 2 
        self.nR = np.array([0, -1, 0])  # normal vector of right face - indexed as 3
        self.nT = np.array([0, 0, 1])  # normal vector of top face - indexed as 4
        self.nBo = np.array([0, 0, -1])  # normal vector of bottom face - indexed as 5
        self.face_normals = np.array([self.nF, self.nB, self.nL, self.nR, self.nT, self.nBo])

        #* faces positions - center points on the planes
        self.rF = np.array([self.box_sizeX / 2, 0, 0])  # Front face position
        self.rB = np.array([-self.box_sizeX / 2, 0, 0]) # Back face position
        self.rL = np.array([0, self.box_sizeY / 2, 0])  # Left face position
        self.rR = np.array([0, -self.box_sizeY / 2, 0])  # Right face position
        self.rT = np.array([0, 0, self.box_sizeZ / 2])  # Top face position
        self.rBo = np.array([0, 0, -self.box_sizeZ / 2])  # Bottom face position
        self.face_centers = np.array([self.rF, self.rB, self.rL, self.rR, self.rT, self.rBo])

    def SphericalToCartesian(self, theta, phi):
        """
        Convert spherical coordinates to Cartesian coordinates.

        Args:
            theta (float): Polar angle in radians.
            phi (float): Azimuthal angle in radians.

        Returns:
            np.ndarray: Cartesian coordinates as a 3-element array.
        """
        x = np.cos(theta)
        y = np.sin(theta) * np.cos(phi)
        z = np.sin(theta) * np.sin(phi)
        return np.array([x, y, z])

    def NewDirectionUnitVector(self, u_old, theta_sc, phi_sc):
        """
        Calculate the new direction unit vector after Compton scattering.
        Note: This part is a bit tricky, because we need to transform the scattering angles from the local coordinate system [where 
        note) the incoming photon is along the x-axis] to the global coordinate system. The transformation is done using a rotation matrix 
        note) based on the original direction of the photon.
        
        Args:
            u_old (np.ndarray): Original direction unit vector.
            theta_sc (float): Scattering polar angle in radians.
            phi_sc (float): Scattering azimuthal angle in radians.

        Returns:
            np.ndarray: New direction unit vector after scattering.
        
        #! Theory behind the transformation:
        * check /simulacije/1-build-up faktor/main/theory/StudenZF/comptonNovaSmer
        . u_old: initial direction unit vector
        . if u_old ~ e_x, then the new direction after scattering is simply given by the spherical to cartesian transformation of the scattering angles (theta_sc, phi_sc)
        . if u_old is not along e_x, we need to rotate the local scattering direction (given by theta_sc, phi_sc) to align with the original direction u_old.
        . the rotation can be done using a rotation matrix R, which is constructed from the original direction u_old and two orthogonal vectors v and w that form a right-handed coordinate system with u_old.
        . the new direction after scattering in the global coordinate system is then given by u_new = R @ u_sc_local, where u_sc_local is the scattering direction in the local coordinate system.        
        """
        
        #u_old = u_old / np.linalg.norm(u_old) # unnecessary as u_old is already a unit vector
        u_sc_local = self.SphericalToCartesian(theta_sc, phi_sc)
        
        # if direction before scattering is close to initial direction of the beam, simple transform
        # --> this is the case of initial incident photon entering the box
        if np.allclose(u_old, [1, 0, 0]):
            return u_sc_local
        
        # if not, we need to do the transformation using rotation matrix
        ref_vec = np.array([1, 0, 0])
        j_local = np.cross(ref_vec, u_old)
        j_local /= np.linalg.norm(j_local) # necessary as ref_vec and u_old are not necessarily orthogonal, so j_local is not necessarily a unit vector!!
        k_local = np.cross(u_old, j_local)
        #k_local /= np.linalg.norm(k_local) # unnecessary as u_old and j_local are already unit vectors and orthogonal, so k_local is also a unit vector!
        
        # local base vectors as columns of the rotation matrix
        B = np.column_stack((u_old, j_local, k_local))
        # apply the rotation to the local scattering direction - this gives us the new direction in the global coordinate system
        u_new = B @ u_sc_local
        # normalize new direction vector
        u_new /= np.linalg.norm(u_new)
        return u_new

    def ComptonScatteringInteraction(self, E, u_particle):
        """
        Simulate a Compton scattering interaction for a photon. Based on the Klein-Nishina formula,
        it samples the scattered photon energy and direction. Draw of new direction is done using 
        rejection sampling method, source: 
        https://geant4.web.cern.ch/documentation/pipelines/master/prm_html/PhysicsReferenceManual/electromagnetic/gamma_incident/compton/compton.html#sampling-the-photon-energy

        Args:
            E (float): Energy of the incoming photon in MeV.
            u_particle (np.ndarray): Direction unit vector of the incoming photon.

        Returns:
            tuple: Scattered photon energy and new direction unit vector.
        """
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
        """
        Calculate the free path length for a photon based on its energy and interaction probabilities.

        Args:
            E (float): Energy of the photon in MeV.

        Returns:
            tuple: Free path length and interaction type ('compton', 'phot', 'pair').
        """
        # Use pre-extracted data instead of accessing dictionary each time
        lac_cs = np.interp(E, self.lac_energy, self.lac_cs)
        lac_pe = np.interp(E, self.lac_energy, self.lac_pe)
        u1, u2 = np.random.uniform(0, 1, 2)
        s_cs = - np.log(u1) / lac_cs
        s_pe = - np.log(u2) / lac_pe
        if E > self.pair_production_threshold:  # Use pre-defined constant
            lac_pp = np.interp(E, self.lac_energy, self.lac_pp)
            u3 = np.random.uniform(0, 1)
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
        """
        Calculate the path length for a particle to exit a 3D box.

        Args:
            r_particle (np.ndarray): Position of the particle.
            u_particle (np.ndarray): Direction unit vector of the particle.

        Returns:
            tuple: Minimum path length to exit the box and the index of the exit plane.

    Bistvo pri računanju: 
    i) Teorija:
        # položaj in smer potovanja delca (unit vector): r0 = [x0,y0,z0] in u = [ux, uy, uz] (oz. omega)
        # enačba ravnine: r * n = d; kjer r=[x,y,z], n=[nx,ny,nz] ~ normalni vektor ploskve, d = skalar = razdalja od izhodišča do ploskve
        # za premik delca velja r(lambda) = r0 + lambda * u; uporabimo enačbo za ravnino: (r0 + lambda * u) * n = d
        # rešitev za lambda: lambda = (d - r0 * n) / (u * n) ~ realno število - razdalja od začetne točke do ploskve, kjer delec seka ploskev
    ii) Implementacija:
        # u * n = cos(alpha); če |u|=|n|=1 & alpha = kot med vektorjema [splošno cos(alpha) = u * n / (|u|*|n|)]  
        # če u * n = 0 --> delec ne seka ploskve, ampak gre vzporedno z njo (smer potovanja točno vzporedna s ploskvijo) --> lambda = inf
        # če u * n > 0 --> delec seka ploskev, ki je v smeri potovanja (npr. če gre proti desni x ploskvi, potem n=[1,0,0] in u*n = ux > 0)
        # če u * n < 0 --> delec seka ploskev, ki je v nasprotni smeri potovanja (npr. če gre proti levi x ploskvi, potem n=[-1,0,0] in u*n = -ux > 0)
        # za vsako ploskev izračunamo lambda, nato pa vzamemo najmanjšo pozitivno lambda, ki nam pove, kje delec zapusti ploskev
        """
        min_path_length = np.inf
        exit_plane_index = -1

        # optional, perhaps more efficient for thick box (num_hvl=d/hvl <~ 2?) - add if statement for the initial direction; we know the 
        # initial path to exit equals box thickness and exit plane is the right one [index 0]
        if np.all(u_particle == [1,0,0]): #could also use np.allclose(u_particle, [1, 0, 0]):
            return self.box_sizeX, 0
        
        for plane_index in range(6):
            #! This line below could be ommited as it is defined for every simulation step - time inefficient, bt more clear code
            r_plane, n_plane = self.face_centers[plane_index], self.face_normals[plane_index]
            #print(f"plane_index:{plane_index}; u_particle:{u_particle}; n_plane={n_plane}")
            relative_direction_scaling_factor = np.dot(u_particle, n_plane)
            #print(f"re_dir_sc_factor:{relative_direction_scaling_factor}")
            #note: if relative_direction_scaling_factor <= 0, it means that the particle is moving away from the plane or parallel to it,
            #note]  so we can skip this plane
            if relative_direction_scaling_factor <= 0:
                continue
            path_length_to_plane = np.abs(np.dot(n_plane, r_particle - r_plane)) / relative_direction_scaling_factor
            if path_length_to_plane < min_path_length:
                min_path_length = path_length_to_plane
                exit_plane_index = plane_index
        return min_path_length, exit_plane_index

    def ParticleStep(self, PreStepParticleInfo):
        """
        Simulate a single step for a particle within the box.

        Args:
            PreStepParticleInfo (dict): Dictionary containing particle information before the step.

        Returns:
            dict: Updated particle information after the step.
        """
        # Extract particle information from the input dictionary
        r_particle = PreStepParticleInfo['r']
        u_particle = PreStepParticleInfo['u']
        E = PreStepParticleInfo['E']


        # Initialize list to store photons to be simulated in the next step
        add_photons_to_simulate = [] 
        # Calculate the path length to exit the box and the interaction free path length
        box_path_length, exit_plane_index = self.BoxPathLength3D(r_particle, u_particle)
        free_path_length, interaction_type = self.PhotonFreePath(E)
        
        #* a) Photon exits the box without interaction
        if free_path_length > box_path_length:
            interaction_type = 'exit'
            E_dep = 0.0
            E_new = E
            u_new = u_particle
            r_new = r_particle + box_path_length * u_particle
            #print(f"   --->Particle exited the box at face index {exit_plane_index} and energy {E_new:.4f} MeV")

        # * b) Photon interacts inside the box
        else:
            r_new = r_particle + free_path_length * u_particle
            exit_plane_index = -1 # set to -1 to indicate interaction inside the box, not exit
            
            if interaction_type == 'compton':
                E_new, u_new = self.ComptonScatteringInteraction(E, u_particle)
                E_dep = E - E_new
                # check for Emin_terminate condition (only necessary for CS)
                if E_new > self.Emin_terminate:
                    # add scattered photon to the list for simulation
                    add_photons_to_simulate.append({
                        'eid': eid,
                        'tid': tid,
                        'E': E_new,
                        'u': u_new,
                        'r': r_new,
                    })
                # if E_new < E_min, then we terminate photon - assumed it is absorbed - need this info for post-processing?
            
            elif interaction_type == 'phot':
                E_new = 0.0
                E_dep = E
                u_new = None
            
            elif interaction_type == 'pair':
                E_new = 0.0
                E_dep = E - self.pair_production_threshold  # Use pre-defined constant
                u_new = None
                # calculate directions of annihilated photons (isotropic distribution)
                phi = np.random.uniform(0, 2 * np.pi)
                th = np.arccos(np.random.uniform(-1, 1))
                x, y, z = np.sin(th) * np.cos(phi), np.sin(th) * np.sin(phi), np.cos(th)
                u_pp = np.array([x, y, z])
                # add the two annihilation photons, emitted in opposite directions, to the list for simulation
                add_photons_to_simulate.extend([
                    {'eid': eid, 'tid': tid + 1, 'E': self.annihilation_energy, 'u': u_pp, 'r': r_new},
                    {'eid': eid, 'tid': tid + 2, 'E': self.annihilation_energy, 'u': -u_pp, 'r': r_new} 
                ])
            
            else:
                raise RuntimeError("Unknown interaction type")
                    
        return {
            'interaction': interaction_type,
            'exit_plane_index': exit_plane_index,
            'E_dep': np.round(E_dep, 6),  # round to 6 decimal places for better readability
            'add_photons_to_simulate': add_photons_to_simulate,
            'E_new': np.round(E_new, 6),  
            'u_new': np.round(u_new, 6) if u_new is not None else None,
            'r_new': np.round(r_new, 6),
                    }

    def simulate_buildup_calc(self, Nsim, report=False):
        """ Simulates the propagation of photons through a defined box geometry.
        Args:
            Nsim (int): Number of simulation events.
            SourcePhotons_template (dict): Template for the source photons, containing:
                - 'E': Energy of the photons (MeV).
                - 'r': Initial position of the photons (for now center left plane: [-a/2,0,0] ~ pencil-beam)
                - 'u': Initial direction unit vector of the photons (for now along the x-axis: [1,0,0]).
        """
        self.Nsim = Nsim
        self.Nout_primary = 0  # number of photons exiting the box w/out interaction
        self.Nout_total = 0  # number of photons exiting the box after interaction(s)
        print(f"Starting simulation: E={self.E0}; num_hvlX = {self.num_hvlX}; {self.Nsim} events\n\n")
        
        # total energy exiting on the right size of the box -> build-up contribution
        self.E_out_tot = 0.0
        # backscattered energy - exit on the face where they entered
        self.E_backscattered = 0.0

        # Printing progress frequency
        #? delete this?
        progress_report_freq = 10  # %
        if Nsim >= 1000:
            eid_progress_set = set(
                int(Nsim * pct / 100) for pct in range(progress_report_freq, 101, progress_report_freq)
            )
        else:
            eid_progress_set = set()

        #* start the clock for simulation time measurement
        start_time = dt.now()

        #* loop over all simulated photons
        for eid in range(Nsim):
            # Print progress
            if eid in eid_progress_set:
                print(f"\tProgress: {eid / Nsim * 100:.2f}%")

            #* add initial photon to the list of photons to simulate - [r, u, E] - position, direction, energy
            photons_to_simulate = [[self.r_entrance,self.u_entrance,self.E0]] # reset the list for new eid

            ##########################################################
            #* loop that continues until photons_to_simulate is empty
            while photons_to_simulate:
                #+ Process the first photon in the list
                # photon info: position, direction, energy
                r_particle, u_particle, E_particle = photons_to_simulate.pop(0)                
                # Calculate the path length to exit the box and the interaction free path length
                box_path_length, exit_plane_index = self.BoxPathLength3D(r_particle, u_particle)
                free_path_length, interaction_type = self.PhotonFreePath(E_particle)
                
                #* a) Photon exits the box without interaction - count Eout; stop tracking
                if free_path_length > box_path_length:
                    if exit_plane_index == 0:  # Assuming exit_plane_index 0 is the right side of the box
                        self.E_out_tot += E_particle
                        self.Nout_total += 1
                        if E_particle == self.E0:  # non-interacted photon
                            self.Nout_primary += 1 
                    elif exit_plane_index == 1:  # Assuming exit_plane_index 1 is the back side of the box
                        self.E_backscattered += E_particle
                    continue  # move to the next photon in the list, if any

                # * b) Photon interacts inside the box
                else:
                    # move particle inside the box to the interaction point
                    r_particle += free_path_length * u_particle
                    
                    #* Compton
                    if interaction_type == 'compton':
                        E_particle, u_particle = self.ComptonScatteringInteraction(E_particle, u_particle)
                        # check for Emin_terminate condition (only necessary for CS)
                        if E_particle > self.Emin_terminate:
                            # add scattered photon to the list for simulation
                            photons_to_simulate.append([r_particle, u_particle, E_particle])
                        #? Do we need this else statement at all here? --> NOPE, if statement is enough
                        #else:
                            # photon energy below termination threshold, stop tracking
                            #continue
                    
                    #* Pair production
                    elif interaction_type == 'pair':
                        # calculate directions of annihilated photons (isotropic distribution)
                        phi = np.random.uniform(0, 2 * np.pi)
                        th = np.arccos(np.random.uniform(-1, 1))
                        x, y, z = np.sin(th) * np.cos(phi), np.sin(th) * np.sin(phi), np.cos(th)
                        u_pp = np.array([x, y, z])
                        # add the two annihilation photons, emitted in opposite directions, to the list for simulation
                        photons_to_simulate.extend([
                            [r_particle, u_pp, self.annihilation_energy],
                            [r_particle, -u_pp, self.annihilation_energy]
                        ])
                    
                    #* Photoeffect 
                    #+ This code is not needed for buildup calculation
                    #elif interaction_type == 'phot':
                     #   continue  # photon is absorbed, so we stop tracking it
                    
                        
                ##########################################################
                ##########################################################

      
        #* end the clock for simulation time measurement
        end_time = dt.now()
        self.simulation_time = end_time - start_time

        # report succesful simulation completion and report results
        print("\nSimulation completed.")
        print(f"Simulation Time:   {self.simulation_time.total_seconds()} sec")
        
        # Code for buildup calculation
        if report:
            print(f"Total energy exiting the box (build-up contribution): {self.E_out_tot:.4f} MeV")
            print(f"Energy of backscattered photons exiting the entrance face: {self.E_backscattered:.4f} MeV")
            print(f"Number of photons exiting the box without interaction: {self.Nout_primary}")
            print(f"Total number of photons exiting the box: {self.Nout_total}")
            # calculate buildup
            E_out_theory = self.Nsim * self.E0 * 2**(-self.num_hvlX)  # theoretical energy exiting the box without buildup
            buildup_factor = self.E_out_tot / E_out_theory
            print(f"Theoretical energy exiting the box without buildup: {E_out_theory:.4f} MeV")
            print(f"Calculated build-up factor: {buildup_factor:.4f}") 

        # OR: return data for postcalculations!
        else:
            return {
                'E_out_tot': np.round(self.E_out_tot, 6),
                'E_backscattered': np.round(self.E_backscattered, 6),
                'Nout_total': self.Nout_total,
                'Nout_primary': self.Nout_primary,
                'simulation_time_sec': self.simulation_time.total_seconds(),
                'E0': self.E0,
                'box_num_hvls': (self.num_hvlX, self.num_hvlY, self.num_hvlZ),  
                'box_size_cm': (self.box_sizeX, self.box_sizeY, self.box_sizeZ),
                'gps_hvl_cm': self.gps_hvl,
                'free_path_length': np.round(free_path_length, 6),
                'box_path_length': np.round(box_path_length, 6),
            }


    