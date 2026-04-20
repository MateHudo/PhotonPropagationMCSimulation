import numpy as np
from datetime import datetime as dt
import matplotlib.pyplot as plt
from collections import deque

class photon_box_propagation_simulator:
    """
    Simulates photon propagation in a 3D box-shaped medium, tracking each step,
    considers photoeffect (PE), Compton scattering (CS), and pair production (PP). Assumes local
    electron energy deposition/absorption and no braking radiation or fluorescence photons (for the latter
    fluorescence yield is important. In low-Z matter like tissue it is ok to neglect it, but in
    high-Z material like Lead, fluorescence yield is high, but those photons are be absorbed locally). 
    Data for attenuation coefficients is gotten from: 
    NIST XCOM (https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html);
    using soft tissue mixture.
    * AIM: to calculate the build-up factor for photons propagating through a box-shaped medium, 
    * and to understand the contribution of different interaction types (PE, CS, PP) to the 
    * energy absorption and photon propagation.
    """
    def __init__(self, box_num_hvls, lac_loader, E0, Emin_terminate=0.001):
        # number of HVLs in x, y, z directions
        self.num_hvlX, self.num_hvlY, self.num_hvlZ = box_num_hvls
        # linear attenuation coefficients loader class instance
        self.lac_loader = lac_loader
        self.lac_energy, self.lac_pp, self.lac_cs, self.lac_pe, self.lac_total = lac_loader.get_lac_data()
        # minimal energy of photons to terminate simulation (unit: MeV)
        self.Emin_terminate = Emin_terminate #! This could be altered and paired with Russ. r. or ommited ...
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
        #print(f"GPS_template: \n\tE0={self.E0}; \n\tr_entr={self.r_entrance}; \n\tu_entr={self.u_entrance}")
        #print("Box HVL (half-value layer) = {:.4f} cm".format(self.half_value_layer))

    def define_box_related_parameters(self):
        self.half_value_layer = np.round(np.log(2) / np.interp(self.E0, self.lac_energy, self.lac_total), 6)  # in cm
        self.box_sizeX = self.num_hvlX * self.half_value_layer  # cm
        self.box_sizeY = self.num_hvlY * self.half_value_layer  # cm
        self.box_sizeZ = self.num_hvlZ * self.half_value_layer  # cm

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

################################################
# HELP METHODS FOR PHOTON PROPAGATION SIMULATION
################################################

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

    def PhotonFreePathForcing(self, E):
        """
        Calculate the free path length for a photon based on its energy and interaction probabilities.
        This method uses variance reduction by forcing interactions, where we force the interaction to always be CS or PP.
        Weighting is applied to account for the bias introduced by forcing the interaction.

        Args:
            E (float): Energy of the photon in MeV.

        Returns:
            tuple: Free path length and interaction type ('compton', 'phot', 'pair').
        """
        # lac_tot
        lac_cs = np.interp(E, self.lac_energy, self.lac_cs)
        lac_pe = np.interp(E, self.lac_energy, self.lac_pe)
        lac_tot = lac_cs + lac_pe
        interaction_type = 'compton'  # default interaction type is CS, but it can be changed to PP if drawn below
        # add pair production if energy is above threshold
        lac_pp = 0
        if E > self.pair_production_threshold:  # Use pre-defined constant
            lac_pp = np.interp(E, self.lac_energy, self.lac_pp)
            lac_tot += lac_pp
            # change int. type to pp if drawn
            if np.random.rand() < lac_pp / (lac_pp + lac_cs):
                interaction_type = 'pair'

        # free path lenght
        free_path_length = - np.log(np.random.rand()) / lac_tot
        # weight
        weight_forcing = (lac_cs + lac_pp) / lac_tot  # weight for forcing CS or PP interaction
        
        return free_path_length, interaction_type, weight_forcing

    def PhotonFreePathCombinedVRM(self, E):
        """
        Calculate the free path length for a photon based on its energy and interaction probabilities.
        This method uses variance reduction methods: pdf manipuation and forcing interactions
        Weighting is applied to account for the bias introduced by forcing the interaction.

        Args:
            E (float): Energy of the photon in MeV.

        Returns:
            tuple: Free path length and interaction type ('compton', 'phot', 'pair').
        """
        # lac_tot
        lac_cs = np.interp(E, self.lac_energy, self.lac_cs)
        lac_pe = np.interp(E, self.lac_energy, self.lac_pe)
        lac_tot = lac_cs + lac_pe
        interaction_type = 'compton'  # default interaction type is CS, but it can be changed to PP if drawn below
        # add pair production if energy is above threshold
        lac_pp = 0
        if E > self.pair_production_threshold:  # Use pre-defined constant
            lac_pp = np.interp(E, self.lac_energy, self.lac_pp)
            lac_tot += lac_pp
            # change int. type to pp if drawn
            if np.random.rand() < lac_pp / (lac_pp + lac_cs):
                interaction_type = 'pair'

        # free path lenght; extend by pef (path extension factor)
        free_path_length = - np.log(np.random.rand()) / lac_tot * self.pef
        # weight forcing
        weight_forcing = (lac_cs + lac_pp) / lac_tot  # weight for forcing CS or PP interaction
        # weight for pdf manipulation - exponential biasing with path extension factor (pef)
        weight_pdf = self.pef * np.exp(-lac_tot * free_path_length * (1 - 1/self.pef))
        # total weight is product of both
        #weight_total = weight_forcing * weight_pdf
        
        return free_path_length, interaction_type, weight_pdf, weight_forcing

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

    def _log(self, message):
        """
        Internal logging function for debugging and tracing the simulation steps.

        Args:
            message (str): The message to be logged.
        
        !! LAHKO BI UVEDEL import logging ...
        """
        if self.verbose:
            print(message)
            #print(f"[DEBUG] {message}")




################################################################################################
# STEP PROCESSING METHODS ETC
################################################################################################

    def process_escaped_photon(self, *args):

        # -------------------
        # UNPACK THE ARGUMENTS
        # -------------------
        exit_plane_index, weight_new, E_particle = args 
        #! alternatives: 
            #* - args -> work with atributes (self; is this better, time consuming??)
            #* -   kwargs - time more consuming as dict? (more readable?)

        # ---------------------------------------------
        # RIGHT SIDE EXIT (non-interacted or interacted)
        # ----------------------------------------------
        if exit_plane_index == 0:  # Assuming exit_plane_index 0 is the right side of the box
            self.E_out_total += weight_new * E_particle
            self.N_out_total += 1
            #self._log(f"Photon exits the box without interaction. E_out = {weight_new * E_particle:.4f} MeV")
            if E_particle == self.E0:  # non-interacted photon
                self.N_out_primaries += 1 
                self.E_out_primaries += weight_new * E_particle
        
        # ------------------
        # BACKSCATTERED EXIT
        # ------------------
        elif exit_plane_index == 1:  # Assuming exit_plane_index 1 is the back side of the box
            self.E_backscattered += weight_new * E_particle
            self.N_backscattered += 1
            #self._log(f"Photon exits the box through backscatter. E_backscattered = {weight_new * E_particle:.4f} MeV")
        
        # ------------
        # LEAKAGE EXIT
        # ------------
        else:
            self.N_leakage += 1
            self.E_leakage += weight_new * E_particle
            #self._log(f"Photon exits the box through a leakage path. E_leakage = {weight_new * E_particle:.4f} MeV (exit plane index: {exit_plane_index})")

    def process_interaction(self, *args):

        # -------------------
        # UNPACK THE ARGUMENTS
        # -------------------
        interaction_type, E_particle, u_particle, photons_to_simulate, r_particle_new = args 
        #note: Velik je teh parametrov, spet dvom a je to smiselno/učinkovito/elegantno...?
        #! alternatives: 
            #* - args -> work with atributes (self; is this better, time consuming??)
            #* -   kwargs - time more consuming as dict? (more readable?)

        #! Need to add returned stuff...
        #? photons_to_sumulate ? 
        #new_photons = []

        # N_absorbed += 1

        # ---------------------------------------------
        # COMPTON SCATTERING
        # ----------------------------------------------
        if interaction_type == 'compton':
            E_particle_new, u_particle_new = self.ComptonScatteringInteraction(E_particle, u_particle)
            self.E_absorbed += E_particle - E_particle_new
            #self._log(f"Compton scattering: E_abs = {(E_particle - E_particle_new):.4f} MeV (E_new = {E_particle_new:.4f} MeV; E_abs_real = {(E_particle - E_particle_new):.4f} MeV)")
            photons_to_simulate.append([r_particle_new, u_particle_new, E_particle_new])
        
        # ---------------- 
        # PAIR PRODUCTION
        # ----------------
        elif interaction_type == 'pair':
            u_pp = self.pair_prodution_new_direction() 
            # add the two annihilation photons, emitted in opposite directions, to the list for simulation
            photons_to_simulate.extend([
                [r_particle_new, u_pp, self.annihilation_energy],
                [r_particle_new, -u_pp, self.annihilation_energy]
            ])
            self.E_absorbed += E_particle - 2 * self.annihilation_energy
            #self._log(f"Pair production: E_dep = {(E_particle - 2 * self.annihilation_energy):.4f} MeV")
        
        # -------------------
        # PHOTOEFFECT (should not occur in forcing method, but we keep it for completeness)
        # -------------------
        else:
            self.E_absorbed += E_particle  # all photon energy is absorbed in the interaction
            #self._log(f"Photoeffect: E_abs = {E_particle:.4f} MeV")
        
        return photons_to_simulate 

    def process_interaction_weighted(self, *args):

        # -------------------
        # UNPACK THE ARGUMENTS
        # -------------------
        interaction_type, E_particle, weight_new, u_particle, photons_to_simulate, r_particle_new = args 
        #note: Velik je teh parametrov, spet dvom a je to smiselno/učinkovito/elegantno...?
        #! alternatives: 
            #* - args -> work with atributes (self; is this better, time consuming??)
            #* -   kwargs - time more consuming as dict? (more readable?)

        #! Need to add returned stuff...
        #? photons_to_sumulate ? 
        #new_photons = []

        # N_absorbed += 1

        # ---------------------------------------------
        # COMPTON SCATTERING
        # ----------------------------------------------
        if interaction_type == 'compton':
            E_particle_new, u_particle_new = self.ComptonScatteringInteraction(E_particle, u_particle)
            self.E_absorbed += weight_new*(E_particle - E_particle_new)  # add absorbed energy to the total
            #self._log(f"Compton scattering: E_abs = {weight_new*(E_particle - E_particle_new):.4f} MeV (E_new = {E_particle_new:.4f} MeV; E_abs_real = {(E_particle - E_particle_new):.4f} MeV)")
            photons_to_simulate.append([r_particle_new, u_particle_new, E_particle_new, weight_new])
        
        # ---------------- 
        # PAIR PRODUCTION
        # ----------------
        elif interaction_type == 'pair':
            u_pp = self.pair_prodution_new_direction() 
            # add the two annihilation photons, emitted in opposite directions, to the list for simulation
            photons_to_simulate.extend([
                [r_particle_new, u_pp, self.annihilation_energy, weight_new],
                [r_particle_new, -u_pp, self.annihilation_energy, weight_new]
            ])
            self.E_absorbed += weight_new*(E_particle - 2 * self.annihilation_energy)
            #self._log(f"Pair production: E_dep = {weight_new*(E_particle - 2 * self.annihilation_energy):.4f} MeV")
        
        # -------------------
        # PHOTOEFFECT (should not occur in forcing method, but we keep it for completeness)
        # -------------------
        #*** LEAVE THIS FOR NOW, problems with E_absorbed, but we dont care for buildup
        # #note: takele zadeve so verjetno problem interpereted jezikov, s C++ bi tole vrjetn elegantno rešu lahko...
        # #note: sreča pa je, da je PE zadnja interakcija vedno, če sploh, ker se foton potem absorbira
        # else:
        #     if self.simulation_method in ['forcing', 'combined_VRM']:
        #         raise ValueError("Photoeffect should not occur in PhotonFreePathForcing, check the code!")
        #     else:
        #         #! For buildup, no need to count E_absorbed actually!
        #         self.E_absorbed += weight_new * E_particle  # all photon energy is absorbed in the interaction
        #         #self._log(f"Photoeffect: E_abs = {weight_new * E_particle:.4f} MeV")
        
        return photons_to_simulate 

    def russian_roulette(self, weight):
        
        # -------------------------------
        # RUSSIAN ROULETTE
        # -------------------------------
        alive = True
        if weight < self.weight_min:
            #self._log(f"\t*Photon weight {weight_new:.4f} below threshold {self.weight_min:.4f}. Applying Russian roulette...")
            if np.random.rand() < self.survival_probability:
                weight /= self.survival_probability
                #self._log(f"Photon survived Russian roulette. New weight: {weight:.4f}")
            else:
                #self._log(f"Photon terminated by Russian roulette.")
                self.N_weight_terminated += 1
                alive = False  # signal that photon is terminated
        
        return alive, weight

    def pair_prodution_new_direction(self):
        phi = np.random.uniform(0, 2 * np.pi)
        th = np.arccos(np.random.uniform(-1, 1))
        return np.array([
            np.sin(th) * np.cos(phi),
            np.sin(th) * np.sin(phi),
            np.cos(th)
            ])

    def apply_path_extension(self,fpl,E):
        """
        
        """
        # extend free path length by pef
        fpl *= self.pef
        # compute weight_pdf for path extension
        lac_tot = np.interp(E, self.lac_energy, self.lac_total)
        weight_pdf = self.pef * np.exp(-lac_tot*fpl*(1-1/self.pef)) # w = p(d)/p_alt(d)
        return fpl, weight_pdf


####################################
# FORCING  FIRST INTERACTION METHODS 
####################################

    def first_interaction_forcing(self):
        """
        Forces the first interaction of a photon inside the box.
        The interaction is CS or PP, we neglect PE as it does not contribute to buildup.
        """
        self.N_steps += 1

        # ----------------------------------------------------------------------------
        # 1. Force the first interaction inside the box - define weight and new position
        # ----------------------------------------------------------------------------
        #todo (test this code and loop succesfull rate, problematic for low N_hvl)
        while True:
            box_path_length, _ = self.BoxPathLength3D(self.r_entrance, self.u_entrance)
            free_path_length, interaction_type, weight_factor = self.PhotonFreePathForcing(self.E0)
            # free path must be less than path to exit
            if free_path_length < box_path_length:
                #self._log(f"\t*Forced interaction: free path length {free_path_length:.4f} cm is less than box path length {box_path_length:.4f} cm.")
                break  # we have a valid forced interaction point, exit the loop
        
        # ---------------------
        # 2. Calculate total weight
        # ---------------------
        weight_total = self.weight_first_interaction * weight_factor
        initial_position = self.r_entrance + free_path_length * self.u_entrance 
        #self._log(f"\t*Free path length drawn: {free_path_length:.4f} cm; new weight: {weight_total:.4f}")

        # ------------------------------------------------------------------------------
        # 3. Define new direction and energy - add photons to the list if int != PE; calculate Edep
        # ------------------------------------------------------------------------------
        return self.process_interaction_weighted(interaction_type, self.E0, weight_total, self.u_entrance, deque([]), initial_position)

    def first_interaction_combined(self):
        """
        Forces the first interaction of a photon inside the box.
        The interaction is CS or PP, we neglect PE as it does not contribute to buildup.
        """
        self.N_steps += 1

        # ----------------------------------------------------------------------------
        # 1. Force the first interaction inside the box - define weight and new position
        # ----------------------------------------------------------------------------
        #todo (test this code and loop succesfull rate, problematic for low N_hvl... or not, contemplte)
        while True:
            #* calculate d_exit and draw fpl, int type
            box_path_length, _ = self.BoxPathLength3D(self.r_entrance, self.u_entrance)
            free_path_length, interaction_type, weight_forcing = self.PhotonFreePathForcing(self.E0)        
            #* prolong path, weight_pdf
            free_path_length, weight_pdf = self.apply_path_extension(free_path_length, self.E0)
            # free path must be less than path to exit
            if free_path_length < box_path_length:
                #self._log(f"\t*Forced interaction: free path length {free_path_length:.4f} cm is less than box path length {box_path_length:.4f} cm.")
                break  
        
        # ---------------------
        # 2. Calculate total weight
        # ---------------------
        #$ weight_force_1st = 1 - 2**(-N_hvl/2) right ??
        weight_total = self.weight_first_interaction * weight_forcing * weight_pdf
        initial_position = self.r_entrance + free_path_length * self.u_entrance  
        #self._log(f"\t*Free path length drawn: {free_path_length:.4f} cm; new weight: {weight_total:.4f}")

        # ------------------------------------------------------------------------------
        # 3. Define new direction and energy - add photons to the list if int != PE; calculate Edep
        # ------------------------------------------------------------------------------
        return self.process_interaction_weighted(interaction_type, self.E0, weight_total, self.u_entrance, deque([]), initial_position)


#########################################
# OLDER SIMULATION CODE (september 2025)
########################################

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
        # fixme: could remove eid and tid for simple simulation of build-up!!
        eid = PreStepParticleInfo['eid']
        tid = PreStepParticleInfo['tid']

        # Initialize list to store photons to be simulated in the next step
        add_photons_to_simulate = [] 
        # Calculate the path length to exit the box and the interaction free path length
        box_path_length, exit_plane_index = self.BoxPathLength3D(r_particle, u_particle)
        free_path_length, interaction_type = self.PhotonFreePath(E)
        
        #* a) Photon exits the box without interaction 
        #fix: it is inefficient to set this condition check as the first one as it is LESS probable (for n_hvl>1...)
        #? How much does this affect t_simulation? Maybe it is insignificant...
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
            #? Is rounding necessaray, time consuming?
            'E_dep': np.round(E_dep, 6),  # round to 6 decimal places for better readability
            'add_photons_to_simulate': add_photons_to_simulate,
            'E_new': np.round(E_new, 6),  
            'u_new': np.round(u_new, 6) if u_new is not None else None,
            'r_new': np.round(r_new, 6),
            
            #* "extras - maybe useful one day // for detailed output or debugging, but not necessary for main simulation results"
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

    def simulate(self, Nsim, doWrite=False):
        """ Simulates the propagation of photons through a defined box geometry.
        Args:
            Nsim (int): Number of simulation events.
        """
        
        self.simulation_method = "simulate"
        self.Nsim = Nsim
        #* counter for photons: they get absorbed or exit the box (on different sides, with different energies, etc.)
        #! - This could be implemented in a more elegant way, for example using a dictionary to store the counts for different categories 
        #!   of photons, but for now let's keep it simple with separate variables for each category!! 
        #! - Could later be ommited and done in post-processing code!! - more efficenent and cleaner simulation code + "decomposition"        
        #! - Can represent a good check of whether the simulation is working correctly: eg. compare number of non-interacted photons exiting
        #!    the box with theoretically expected (N_out_non_interacted/Nsim = 2**(-num_hvlX)) - this is necessary, but NOT SUFFICIENT condition ofcourse
        self.N_steps = 0  # number of steps simulated (for all photons = number of lines in the output file= number of ParticleStep calls)
        self.N_out_primaries = 0  # number of photons exiting the box ON RIGHT SIDE w/out interaction
        self.N_out_total = 0  # total number of photons exiting the box ON RIGHT SIDE (both interacted and non-interacted)
        self.N_backscattered = 0  # number of photons backscattered (exit on the face where they entered)
        self.N_leakage = 0  # number of photons exiting on the four side faces of the box
        # same for energies
        self.E_out_total = 0.0 # right side
        self.E_absorbed = 0.0
        self.E_backscattered = 0.0
        self.E_leakage = 0.0


        print(f"Starting simulation: E={self.E0} MeV; num_hvlX = {self.num_hvlX}; {self.Nsim} events\n\n")
        
        # create a timestamped filename for saving simulation info, store inside "initial_output" folder
        timestamp_str = dt.now().strftime("%Y%m%d_%H%M%S")
        output_file= f"initial_output/run_{timestamp_str}_{self.E0}MeV_{self.num_hvlX}HVL.txt"

        #define a writer function depending on doWrite
        #? Is this the best way to do it? maybe it is better to write the file in batches, for example after every 1000 events, 
        #? to reduce the number of file writes and improve performance? [but for now, let's keep it simple and write line by line, as it is 
        #? easier to implement and understand]
        if doWrite:
            def write_line(msg):
                with open(output_file,"a") as f:
                    f.write(msg)
        else:
            def write_line(msg):
                pass # no-op

        # Printing progress frequency
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
                print(f"\tProgress: {eid / Nsim * 100:.2f}% ({(dt.now() - start_time).total_seconds():.2f} sec)")

            # add initial photon to simulate - if eid+tid is ommited, this cn be defined once instead of for every eid...
            photons_to_simulate = [{
                'eid': eid,
                'tid': 0,
                'r': self.r_entrance,
                'u': self.u_entrance,
                'E': self.E0,
            }]
    

            #* loop that continues until photons_to_simulate is empty
            photons_to_simulate = deque([{
                'eid': eid,
                'tid': 0,
                'r': self.r_entrance,
                'u': self.u_entrance,
                'E': self.E0,
            }])
    

            #* loop that continues until photons_to_simulate is empty
            while photons_to_simulate:
                # Process the first photon in the list
                PreStepData = photons_to_simulate.popleft()
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
                

                #* write the result to the file
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


                ##* Build-up contribution, backscatter, leakage, absorbed energy, ...
                #! This could also be moved to post-processing code... [simulation time etc] - or ommited at all!!
                # First counting all simulation steps
                self.N_steps += 1
                # a) energy absorbed inside the box
                if exit_plane_index == -1:
                    self.E_absorbed += Edep
                # b) build-up contribution
                elif exit_plane_index == 0:  # Assuming exit_plane_index 0 is the right side of the box
                    self.N_out_total += 1
                    self.E_out_total += Enew
                    if Enew == self.E0:
                        self.N_out_primaries += 1
                # c) backscattered photons - exit on the face where they entered
                elif exit_plane_index == 1:  # Assuming exit_plane_index 1 is the back side of the box
                    self.N_backscattered += 1
                    self.E_backscattered += Enew
                # d) leakage photons - exit on the side faces of the box
                else:
                    self.N_leakage += 1
                    self.E_leakage += Enew
        
        #* end the clock for simulation time measurement
        end_time = dt.now()
        self.simulation_time = end_time - start_time

        # report succesful simulation completion and report results
        print("\nSimulation completed.")
        print(f"Simulation Time:   {self.simulation_time.total_seconds()} sec")
        # calculate buildup (this could be erased once i use code for different purposes than buildup calculation)
        E_out_theory = self.Nsim * self.E0 * 2**(-self.num_hvlX)
        buildup = self.E_out_total / E_out_theory
        print(f"Buildup Factor:    {buildup:.6f} (E_out={self.E_out_total:.6f} MeV; E_out_theory={E_out_theory:.6f} MeV)")
        #* final report of the results
        #self.report_results()
 


################################################################################################
# SIMULATION METHODS: normal (buildup), pfd manipulation, forcing, combined VRM
################################################################################################

    def simulate_buildup_calc(self):
        """ Simulates the propagation of photons through a defined box geometry.
            #! This is a simplified version of the simulate() method, focused only on calculating the buildup 
            #! contribution (energy exiting on the right side of the box). Most important change seems to be actually putting
            #! PhotonStep() and simulate() methods together - aim (achieved, seemingly by ~30%) was faster simulation time. 
            #! Importantly: i get same results so code seems o be correct!
            #! All further anaylsis - comparing codes etc - is for later
        Args:
            Nsim (int): Number of simulation events.
        """

        # ----------
        # ...
        # ----------
        for eid in range(self.Nsim):
            #self._log(f"\nSimulating event {eid+1}/{self.Nsim}...")
            # Print progress
            #if eid in self.eid_progress_set:
            #    print(f"\tProgress: {eid / self.Nsim * 100:.2f}% ({(dt.now() - self.start_time).total_seconds():.2f} sec)")

            #* add initial photon to the list of photons to simulate - [r, u, E] - position, direction, energy
            photons_to_simulate = deque([[self.r_entrance,self.u_entrance,self.E0]]) # reset the list for new eid

            # ---------------
            # loop while photons in the list to simulate
            # ---------------
            while photons_to_simulate:
                #+ Process the first photon in the list
                self.N_steps += 1
                # photon info: position, direction, energy
                r_particle, u_particle, E_particle = photons_to_simulate.popleft()     
                #self._log(f"\t*Processing photon with E={E_particle:.4f} MeV at position {r_particle} and direction {u_particle}...")           
                # Calculate the path length to exit the box and the interaction free path length
                box_path_length, exit_plane_index = self.BoxPathLength3D(r_particle, u_particle)
                free_path_length, interaction_type = self.PhotonFreePath(E_particle)
                
                # --------------------
                # ANALYZE STEP OUTCOME
                # --------------------
                #* a) Photon interacts inside the box (more likely, so we check this first)
                if free_path_length < box_path_length:
                    # move particle inside the box to the interaction point
                    r_particle_new = r_particle + free_path_length * u_particle
                    self.N_box += 1
                    photons_to_simulate = self.process_interaction(interaction_type, E_particle, u_particle, photons_to_simulate, r_particle_new)
                #* b) Photon exits the box without interaction - count Eout; use old weight, stop tracking
                else:
                    self.process_escaped_photon(exit_plane_index, 1, E_particle)

    def simulate_pdf_manipulation(self):
        """
        Simulates the propagation ...
        """

        # ---------------
        # ...
        # ---------------
        for eid in range(self.Nsim):
            #self._log(f"\nSimulating event {eid+1}/{self.Nsim}")
            # Print progress
            if eid in self.eid_progress_set:
                print(f"\tProgress: {eid / self.Nsim * 100:.2f}% ({(dt.now() - self.start_time).total_seconds():.2f} sec)")

            # add initial photon
            photons_to_simulate = deque([[self.r_entrance,self.u_entrance,self.E0, self.start_weight]]) # reset the list for new eid

            # ---------------------
            # ...
            # ---------------------
            while photons_to_simulate:
                # ------------------------------------
                # Process the first photon in the list
                # ------------------------------------
                self.N_steps += 1
                r_particle, u_particle, E_particle, weight = photons_to_simulate.popleft()     
                #self._log(f"\t*Processing photon with E={E_particle:.4f} MeV at position {r_particle}, direction {u_particle} and weight {weight:.4f}")           
                box_path_length, exit_plane_index = self.BoxPathLength3D(r_particle, u_particle)
                
                # --------------------------------
                # Free path (extended); new weight
                # --------------------------------
                free_path_length, interaction_type = self.PhotonFreePath(E_particle)
                free_path_length, weight_pdf = self.apply_path_extension(free_path_length, E_particle)
                #note: use energy-dependent HVL (not HVL at E0) so secondary photons (E != E0) get the correct step weight
                weight_new = weight * weight_pdf 
                #self._log(f"\t*Free path lenght drawn: {free_path_length:.4f} cm; new weight: {weight_new:.4f}")
                
                # -------------------------------
                # RUSSIAN ROULETTE
                # -------------------------------
                alive, weight_new = self.russian_roulette(weight_new)
                #* terminate photon if it did not survive Russian roulette
                if not alive:
                    continue
                
                # --------------------
                # ANALYZE STEP OUTCOME
                # --------------------
                #* a) Photon interacts inside the box (more likely, so we check this first)
                if free_path_length < box_path_length:
                    # move particle inside the box to the interaction point
                    r_particle_new = r_particle + free_path_length * u_particle
                    self.N_box += 1
                    photons_to_simulate = self.process_interaction_weighted(interaction_type, E_particle, weight_new, u_particle, photons_to_simulate, r_particle_new)
                #* b) Photon exits the box without interaction - count Eout; use old weight, stop tracking
                else:
                    self.process_escaped_photon(exit_plane_index, weight_new, E_particle)

    def simulate_forcing_interactions(self):
        """
        Comments...
        """

        # -------------------------------
        # LOOP OVER ALL SIMULATED PHOTONS
        # --------------------------------
        for eid in range(self.Nsim):
            #self._log(f"\nSimulating event {eid+1}/{self.Nsim}")
            # Print progress
            if eid in self.eid_progress_set:
                print(f"\tProgress: {eid / self.Nsim * 100:.2f}% ({(dt.now() - self.start_time).total_seconds():.2f} sec)")

            # --------------------
            # INITIALIZE NEW PHOTON
            # --------------------
            if self.force_first_interaction:
                #self._log(f"\t*Processing the initial photon with E={self.E0:.4f} MeV at position {self.r_entrance} and direction {self.u_entrance} with forcing the first interaction...")
                photons_to_simulate = self.first_interaction_forcing()
            else:
                photons_to_simulate = deque([[self.r_entrance,self.u_entrance,self.E0, self.start_weight]]) # reset the list for new eid

            # -----------------------------
            # LOOP OVER PHOTONS TO SIMULATE
            # -----------------------------
            while photons_to_simulate:

                # ------------------------------------
                # PROCESS THE FIRST PHOTON IN THE LIST
                # ------------------------------------
                self.N_steps += 1
                r_particle, u_particle, E_particle, weight = photons_to_simulate.popleft() 
                #self._log(f"\t*Processing photon with E={E_particle:.4f} MeV at position {r_particle}, direction {u_particle} and weight {weight:.4f}")           
                box_path_length, exit_plane_index = self.BoxPathLength3D(r_particle, u_particle)
                free_path_length, interaction_type, weight_factor = self.PhotonFreePathForcing(E_particle)
                
                # --------------------
                # ANALYZE STEP OUTCOME
                # --------------------
                if free_path_length < box_path_length:
                    # move particle inside the box to the interaction point
                    r_particle_new = r_particle + free_path_length * u_particle
                    self.N_box += 1

                    # -------------------------------
                    # UPDATE WEIGHT, RUSSIAN ROULETTE
                    # -------------------------------
                    weight_new = weight * weight_factor
                    #self._log(f"\t*Free path length drawn: {free_path_length:.4f} cm; new weight: {weight_new:.4f}")
                    #* Russian Roulette
                    alive, weight_new = self.russian_roulette(weight_new)
                    if not alive:
                        continue
                    
                    # ----------------
                    # PROCESS INTERACTION
                    # ----------------
                    photons_to_simulate = self.process_interaction_weighted(interaction_type, E_particle, weight_new, u_particle, photons_to_simulate, r_particle_new)
                
                # -------------
                # EXIT FROM BOX
                #--------------
                else:
                    self.process_escaped_photon(exit_plane_index, weight, E_particle) #! Here must use weight, not weight_new

    def simulate_combinedVRM(self):
        """
        Comments...
        """

        # -------------------------------
        # LOOP OVER ALL SIMULATED PHOTONS
        # -------------------------------
        for eid in range(self.Nsim):

            # -------------------------
            # #### PRINT PROGRESS #####
            # -------------------------
            #self._log(f"\nSimulating event {eid+1}/{self.Nsim}")
            if eid in self.eid_progress_set:
                print(f"\tProgress: {eid / self.Nsim * 100:.2f}% ({(dt.now() - self.start_time).total_seconds():.2f} sec)")

            # --------------------
            # INITIALIZE NEW PHOTON
            # --------------------
            if self.force_first_interaction:
                #self._log(f"\t*Processing the initial photon with E={self.E0:.4f} MeV at position {self.r_entrance} and direction {self.u_entrance} with forcing the first interaction and path length manipulation...")
                photons_to_simulate = self.first_interaction_combined()
            else:
                photons_to_simulate = deque([[self.r_entrance,self.u_entrance,self.E0, self.start_weight]])


            # -----------------------------
            # LOOP OVER PHOTONS TO SIMULATE
            # -----------------------------
            while photons_to_simulate:

                # ------------------------------------
                # PROCESS THE FIRST PHOTON IN THE LIST
                # ------------------------------------
                self.N_steps += 1
                r_particle, u_particle, E_particle, weight = photons_to_simulate.popleft() 
                #self._log(f"\t*Processing photon with E={E_particle:.4f} MeV at position {r_particle}, direction {u_particle} and weight {weight:.4f}")           
                box_path_length, exit_plane_index = self.BoxPathLength3D(r_particle, u_particle)
                free_path_length, interaction_type, weight_pdf, weight_forcing = self.PhotonFreePathCombinedVRM(E_particle)
                
                
                # --------------------
                # ANALYZE STEP OUTCOME
                # --------------------
                #* a) Photon interacts inside the box (more likely, so we check this first)
                if free_path_length < box_path_length:
                    # move particle inside the box to the interaction point
                    r_particle_new = r_particle + free_path_length * u_particle
                    self.N_box += 1
                    # w_new and Russian Roulette
                    weight_new = weight * weight_pdf * weight_forcing
                    alive, weight_new = self.russian_roulette(weight_new)
                    if not alive:
                        continue
                    photons_to_simulate = self.process_interaction_weighted(interaction_type, E_particle, weight_new, u_particle, photons_to_simulate, r_particle_new)
                #* b) Photon exits the box without interaction - count Eout; use old weight, stop tracking
                else:
                    # w_new (!! without weight_forcing !!) and RR
                    weight_new = weight * weight_pdf
                    alive, weight_new = self.russian_roulette(weight_new)
                    if not alive:
                        continue
                    self.process_escaped_photon(exit_plane_index, weight_new, E_particle)





#####################################################
# CONFIGURATION (SETUP)
#####################################################

    def read_config(self,config):
        """
        Read the simulation configuration and set up the necessary parameters.

        Args:
            config (dict): Dictionary containing simulation configuration parameters.
        """
        self.config = config
        #self._log(f"Configuration read: {config}") 

        self.simulation_method = config['simulation_method']
        self.Nsim = config['Nsim']
        self.Srep_index = config['Srep_index'] # simulation repetition index
        self.verbose = config['verbose']
        self.start_weight = config['start_weight']
        self.weight_min = config['weight_min'] # threshold for Russian roulette (to terminate low-weight photons and avoid infinite tracking)
        self.survival_probability = config['survival_probability']  # survival probability for Russian roulette (if photon weight falls below threshold, it has this probability to survive and continue tracking with boosted weight)
        self.pef = config['path_extension_factor']  # factor by which the free path
        self.force_first_interaction = config['force_first_interaction']  # whether to force the first interaction of each photon (for variance reduction)

        #* place for additional attributes definitios
        self.weight_first_interaction = 1 - 2**(-self.num_hvlX)
        if self.simulation_method == "combined":
            self.weight_first_interaction = 1 - 2**(-self.num_hvlX/self.pef)
        #? kako bi tole bolš naredu? Nekateri parameteri so relavantni le pri nekaterih metodah... (weight_min, pef, ...)
        self.simulate_primaries = True
        if self.simulation_method in ["forcing", "combined"] and self.force_first_interaction:
            self.simulate_primaries = False


#################
#################
# RUN SIMULATION
#################
#################

    def run(self, config):
        """
        Run the simulation based on the provided configuration.
        """
        # -------------------------------------
        # READ CONFIG, RESET COUNTING VARIABLES
        # -------------------------------------
        self.read_config(config) 
        self.reset_counting()

        # ---------------------------------------------------------------------------
        # PRINT SIMULATION SETTINGS; start the clock for simulation time measurement
        # ---------------------------------------------------------------------------
        #fix: This shall maybe also be managed, changed and/or packed in its own method...
        #print(f"\n\nStarting simulation.")
        #print(f"Setup: method={self.simulation_method}, E0={self.E0} MeV, num_hvlX = {self.num_hvlX}, {self.Nsim} events")
        #print(f"\t?Force first interaction: {self.force_first_interaction}") ##? what to do with this?? 
        

        # ----------------------------
        # START SIMULATION TIME MEASUREMENT
        # ----------------------------
        self.start_time = dt.now()

        # --------------------------------------------
        # RUN SIMULATION (based on the selected method)
        # --------------------------------------------
        #note: can use different mathod names for some methods! Suggestion: all methods manes of same length - more elegant results!
        if self.simulation_method in ['combined','combine',"combined_vrm"]:
            self.simulate_combinedVRM()
        elif self.simulation_method == 'forcing':
            self.simulate_forcing_interactions()
        elif self.simulation_method in ['pdf', 'pdf_man']:
            self.simulate_pdf_manipulation()
        elif self.simulation_method == 'buildup':
            self.simulate_buildup_calc()
        else:
            raise ValueError(f"Unknown simulation method: {self.simulation_method}")

        
        # -------------------------------
        # END SIMULATION TIME MEASUREMENT
        # -------------------------------
        self.end_time = dt.now()
        self.simulation_time = self.end_time - self.start_time

        # -------------------------------------
        # SIMULATION COMPLETED - REPORT RESULTS
        # -------------------------------------
        #print("\nSimulation completed.")
        #print(f"Simulation Time:   {self.simulation_time.total_seconds()} sec")
        
        # calculate buildup
        buildup_factor = self.E_out_total / (self.Nsim * self.E0 * 2**(-self.num_hvlX))
        if not self.simulate_primaries:
            buildup_factor += 1 # add primaries contribution
        #print(f"Calculated build-up factor: {buildup_factor:.4f}") 
        
        # post simulation counting
        self.post_simulation_counting(report=False)
        # simulation info
        return self.simulation_info()



#################################
# COUNTING PHOTONS, REPORT AND SAVE
#################################

    def reset_counting(self):
        """
        Comments...
        """

        #* counting photons and energies (shall focus on buildup clacultion + checking if code works properly - compare to simulate method)
        self.N_steps = 0  # photon step is defined as one photon discrete move (or attempted move)
        self.N_steps_successful = 0  #$ N_steps - N_terminated_total - a bit non-necessary...
        
        # out
        #? better name: transmitted??
        #$ "total" = "primaries" + "secondaries"
        self.N_out_total = 0  # total number of photons exiting the box ON RIGHT SIDE (both interacted and non-interacted)
        self.N_out_primaries = 0 # non-interacted
        self.N_out_secondaries = 0 # those that intearct + annihilation photons from pair production - buildup contribution
        self.E_out_total = 0.0
        self.E_out_primaries = 0.0
        self.E_out_secondaries = 0.0
        # backscatter
        self.N_backscattered = 0  # number of photons backscattered (exit on the face where they entered)
        self.E_backscattered = 0.0
        # leakage (side faces)
        self.N_leakage = 0  # number of photons exiting on the four side faces of the box
        self.E_leakage = 0.0
        # total exit
        self.N_exit = 0
        self.E_exit = 0.0
        
        #! OD TUKI NAPREJ DOL SO SAMI DODATKI, "JAJCA" PO DOMAČ, PUST TO ZDEJ...
        # box - absorbed by the medium or termination based on energy or weight threshold
        #$ "box" = "absorbed" + "terminated by weight threshold" + "terminated by energy threshold"
        self.N_box = 0  # number of photons absorbed in the medium (photoeffect or termination, based on w_min or E_min threshold)
        self.E_box = 0.0
        self.N_absorbed = 0  # number of photons absorbed in the medium (photoeffect only)
        self.E_absorbed = 0.0
        self.N_weight_terminated = 0  # number of photons terminated by Russian roulette
        self.N_energy_terminated = 0  # number of photons terminated by energy threshold
        self.N_terminated_total = 0
        #? those three below needed???
        self.E_weight_terminated = 0.0  # energy of photons terminated by Russian roulette
        self.E_energy_terminated = 0.0  # energy of photons terminated by energy threshold
        self.E_terminated_total = 0.0
        # counting interactions
        self.N_interactions_total = 0
        self.N_compton = 0
        self.N_photoeffect = 0
        self.N_pair_production = 0
        # buildup factor - initialize, not known until the end of the simulation when we have all the counting results
        self.buildup_factor = 0.0

        # Printing progress frequency
        self.progress_report_freq = self.config["progress_report_frequency"]  # %
        if self.Nsim >= 1000:
            self.eid_progress_set = set(
                int(self.Nsim * pct / 100) for pct in range(self.progress_report_freq, 101, self.progress_report_freq)
            )
        else:
            self.eid_progress_set = set()

    def post_simulation_counting(self, report=False):
        """
        Report the results of the counting after the simulation is complete.
        """
        if not hasattr(self, 'N_steps'):
            raise ValueError("Calling post_simulation_counting denied. Please run the simulation before calling post_simulation_counting().")

        #* Calculate non-counted photons
        # photons out
        self.N_out_secondaries = self.N_out_total - self.N_out_primaries
        self.E_out_secondaries = self.E_out_total - self.E_out_primaries
        # photons exit = out + leak + backscatter
        self.N_exit = self.N_out_total + self.N_leakage + self.N_backscattered
        self.E_exit = self.E_out_total + self.E_leakage + self.E_backscattered
        # termination and absorption
        self.N_terminated_total = self.N_weight_terminated + self.N_energy_terminated 
        self.E_terminated_total = self.E_weight_terminated + self.E_energy_terminated 
        self.N_steps_successful = self.N_steps - self.N_terminated_total
        self.N_absorbed = self.N_box - self.N_terminated_total
        self.E_absorbed = self.N_absorbed * self.E0  # assuming absorbed ph. deposit their full energy - wrong for VRM (weights)
        self.E_box = self.E_absorbed + self.E_terminated_total  # total energy absorbed in the box (absorbed + terminated)
        # interactions
        self.N_interactions_total = self.N_steps - self.N_exit - self.N_terminated_total
        # cs, pp, pe...

        # counted vs theortical expected
        self.N_processed_theory = self.Nsim + self.N_pair_production #! wrong if E>E_pp_threshold and we dont count pp events
        self.N_processed_simulation = self.N_exit + self.N_box
        # self.N_processed_difference = self.N_processed_theory - self.N_processed_simulation
        self.E_total_simulation = self.Nsim * self.E0
        self.E_processed_simulation = self.E_exit + self.E_box # need to cound E_box
        # self.E_simulation_diff = self.E_total_simulation - self.E_processed_simulation

        # buildup
        self.N_out_primaries_theory = self.Nsim * 2**(-self.num_hvlX)  # expected number of primaries out
        self.E_out_primaries_theory = self.N_out_primaries_theory * self.E0  # expected energy of primaries out
        self.buildup_factor = self.E_out_total / self.E_out_primaries_theory  # B = E_out_total / E_out_theory
        # add 1 if 1st interacion is forced and we don't simulate primaries
        if not self.simulate_primaries:
            self.buildup_factor += 1

        if report:
            # Print formatted results
            print("\n" + "="*70)
            print(f"{'SIMULATION COUNTING RESULTS':^70}")
            print("="*70)

            print(f"\n{'SIMULATION / EXECUTION':^70}")
            print("-"*70)
            print(f"  Simulation method:            {self.simulation_method}")
            print(f"  Incident energy E0:           {self.E0:.6f} MeV")
            print(f"  Number of primary photons:    {self.Nsim}")
            print(f"  Shield thickness:             {self.num_hvlX} HVL")
            if hasattr(self, 'simulation_time'):
                print(f"  Simulation time:              {self.simulation_time.total_seconds():.2f} sec")

            print(f"\n{'PHOTON COUNTS':^70}")
            print("-"*70)
            print(f"  Total simulation steps:       {self.N_steps}")
            print(f"  Successful steps:             {self.N_steps_successful}")
            print(f"  Total interactions:           {self.N_interactions_total}")
            print(f"  Compton interactions:         {self.N_compton}")
            print(f"  Photoeffect interactions:     {self.N_photoeffect}")
            print(f"  Pair production interactions: {self.N_pair_production}")

            print(f"\n{'EXITED PHOTONS':^70}")
            print("-"*70)
            print(f"  Out total:                    {self.N_out_total}")
            print(f"  Out primaries:                {self.N_out_primaries}")
            print(f"  Out secondaries:              {self.N_out_secondaries}")
            print(f"  Backscattered:                {self.N_backscattered}")
            print(f"  Leakage:                      {self.N_leakage}")
            print(f"  Total exited photons:         {self.N_exit}")

            print(f"\n{'PHOTONS ENDING IN SHIELD':^70}")
            print("-"*70)
            print(f"  In box total:                 {self.N_box}")
            print(f"  Absorbed:                     {self.N_absorbed}")
            print(f"  Weight terminated:            {self.N_weight_terminated}")
            print(f"  Energy terminated:            {self.N_energy_terminated}")
            print(f"  Total terminated:             {self.N_terminated_total}")

            print(f"\n{'ENERGY ACCOUNTING':^70}")
            print("-"*70)
            print(f"  E_out total:                  {self.E_out_total:.6f} MeV")
            print(f"  E_out primaries:              {self.E_out_primaries:.6f} MeV")
            print(f"  E_out secondaries:            {self.E_out_secondaries:.6f} MeV")
            print(f"  E_backscattered:              {self.E_backscattered:.6f} MeV")
            print(f"  E_leakage:                    {self.E_leakage:.6f} MeV")
            print(f"  E_exit total:                 {self.E_exit:.6f} MeV")
            print(f"  E_absorbed:                   {self.E_absorbed:.6f} MeV")
            print(f"  E_weight_terminated:          {self.E_weight_terminated:.6f} MeV")
            print(f"  E_energy_terminated:          {self.E_energy_terminated:.6f} MeV")
            print(f"  E_terminated total:           {self.E_terminated_total:.6f} MeV")
            print(f"  E_box total:                  {self.E_box:.6f} MeV")

            print(f"\n{'BALANCE CHECKS':^70}")
            print("-"*70)
            print(f"  N_processed theory:           {self.N_processed_theory}")
            print(f"  N_processed simulation:       {self.N_processed_simulation}")
            print(f"  N difference:                 {self.N_processed_theory - self.N_processed_simulation}")
            print(f"  E_total simulation:           {self.E_total_simulation:.6f} MeV")
            print(f"  E_processed simulation:       {self.E_processed_simulation:.6f} MeV")
            print(f"  E difference:                 {self.E_total_simulation - self.E_processed_simulation:.6f} MeV")

            print(f"\n{'BUILD-UP FACTOR':^70}")
            print("-"*70)
            print(f"  E_out theory:                 {self.Nsim * self.E0 * 2**(-self.num_hvlX):.6f} MeV")
            print(f"  B = E_out/E_out_theory:       {self.buildup_factor:.6f}")
            print("="*70)

    def simulation_info(self):
        """
        Comments...
        """
        if not hasattr(self, 'Nsim'):
            raise RuntimeError("Simulation not yet run. Please run simulate() before getting results.")

        return {
            # most important
            "method": self.simulation_method,
            "E0": self.E0,
            "n_hvl_x": self.num_hvlX,
            "Srep_index": self.Srep_index,
            "simulation_time_sec": self.simulation_time.total_seconds(),
            "N_stps_per_Nsim": self.N_steps / self.Nsim,
            "buildup_factor": self.buildup_factor,

            
            # other simulation parameters
            "n_hvl_y": self.num_hvlY,
            "n_hvl_z": self.num_hvlZ,
            "Nsim": self.Nsim,
            
            # additional, not used in every method
            "force_first_interaction": self.force_first_interaction if self.simulation_method in ["forcing", "combined"] else None,
            "path_extension_factor": self.pef if self.simulation_method in ["pdf", "combined"] else None,
            "weight_min": self.weight_min if self.simulation_method in ["forcing", "combined", "pdf"] else None,
            "survival_probability": self.survival_probability if self.simulation_method in ["forcing", "combined", "pdf"] else None,

            
            # photon fate distribution
            "N_out_total": self.N_out_total,
            "E_out_total": self.E_out_total,
            "N_out_primaries": self.N_out_primaries,
            "E_out_primaries": self.E_out_primaries,
            "N_out_secondaries": self.N_out_secondaries,
            "E_out_secondaries": self.E_out_secondaries,
            # additionally...
            "N_backscattered": self.N_backscattered,
            "E_backscattered": self.E_backscattered,
            "N_leakage": self.N_leakage,
            "E_leakage": self.E_leakage,
            "N_absorbed": self.N_absorbed,
            "E_absorbed": self.E_absorbed,

            # calculated parameters
            "N_steps": self.N_steps,
            "N_stps_per_sec": self.N_steps / self.simulation_time.total_seconds(),
            "Nsim_per_sec": self.Nsim / self.simulation_time.total_seconds(),

            # other...
            #
        }



##############################
# Older report method ...
    def report_and_save(self, save_to_file=False, print_report=True):
        """
        Build and optionally print/save formatted simulation results.

        Args:
            save_to_file (bool): If True, save report text to a file.
            print_report (bool): If True, print report to console.

        Returns:
            str: Formatted report text.
        """
        if not hasattr(self, 'Nsim'):
            raise RuntimeError("Simulation not yet run. Please run simulate() before reporting results.")

        # --------------------------------------------
        #  SOME ADDITIONAL CALCULATIONS FOR THE REPORT
        # --------------------------------------------
        #$ look at report_counting_results() method
        
        # calculate steps executed per second

        #* Build formatted simulation report text
        report_lines = []
        report_lines.append("")
        report_lines.append("=" * 70)
        report_lines.append(f"{'PHOTON PROPAGATION SIMULATION REPORT':^70}")
        report_lines.append("=" * 70)

        report_lines.append(f"\n{'SIMULATION PARAMETERS':^70}")
        report_lines.append("-" * 70)
        report_lines.append(f"  Incident Energy (E0):       {self.E0:.4f} MeV")
        report_lines.append(f"  Box Dimensions (HVLs):      X={self.num_hvlX}, Y={self.num_hvlY}, Z={self.num_hvlZ}")
        report_lines.append(f"  Box Dimensions (cm):        X={self.box_sizeX:.4f}, Y={self.box_sizeY:.4f}, Z={self.box_sizeZ:.4f}")
        report_lines.append(f"  Method:                     {self.simulation_method}")
        report_lines.append(f"  Number of Events:           {self.Nsim}")
        report_lines.append(f"  Total Simulated Energy:     {self.E_total_simulation:.6f} MeV")

        report_lines.append(f"\n{'SIMULATION EXECUTION':^70}")
        report_lines.append("-" * 70)
        report_lines.append(f"  Total Simulation Time:      {self.simulation_time.total_seconds():.2f} sec ({self.simulation_time.total_seconds()/60:.2f} min)")
        report_lines.append(f"  Total Simulation Steps:     {self.N_steps}")
        report_lines.append(f"  Average Steps per Event:    {self.N_steps/self.Nsim:.2f}")
        report_lines.append(f"  Average Events per Time:    {self.Nsim / self.simulation_time.total_seconds():.2f} events/sec")
        report_lines.append(f"  Simulation Efficiency:      {self.N_steps / self.simulation_time.total_seconds():.2f} steps/sec")

        report_lines.append(f"\n{'PHOTON FATE DISTRIBUTION':^70}")
        report_lines.append("-" * 80)
        report_lines.append(f"  {'Category':<30} {'Count':<8} {'%':<10} {'Energy (MeV)':<15} {'%':<10}")
        report_lines.append("-" * 80)
        
        # Transmitted (exiting right side)
        pct_out_count = (self.N_out_total / self.Nsim * 100) if self.Nsim > 0 else 0
        pct_out_energy = (self.E_out_total / self.E_total_simulation * 100) if self.E_total_simulation > 0 else 0
        report_lines.append(f"  {'Transmitted (Right)':<30} {self.N_out_total:<8} {pct_out_count:<10.2f} {self.E_out_total:<15.6f} {pct_out_energy:<10.2f}")
        
        # Backscattered (exiting left side)
        pct_back_count = (self.N_backscattered / self.Nsim * 100) if self.Nsim > 0 else 0
        pct_back_energy = (self.E_backscattered / self.E_total_simulation * 100) if self.E_total_simulation > 0 else 0
        report_lines.append(f"  {'Backscattered (Left)':<30} {self.N_backscattered:<8} {pct_back_count:<10.2f} {self.E_backscattered:<15.6f} {pct_back_energy:<10.2f}")
        
        # Leakage (exiting through sides)
        pct_leak_count = (self.N_leakage / self.Nsim * 100) if self.Nsim > 0 else 0
        pct_leak_energy = (self.E_leakage / self.E_total_simulation * 100) if self.E_total_simulation > 0 else 0
        report_lines.append(f"  {'Leakage (Sides)':<30} {self.N_leakage:<8} {pct_leak_count:<10.2f} {self.E_leakage:<15.6f} {pct_leak_energy:<10.2f}")
        
        # Absorbed
        pct_abs_count = (self.N_absorbed / self.Nsim * 100) if self.Nsim > 0 else 0
        pct_abs_energy = (self.E_absorbed / self.E_total_simulation * 100) if self.E_total_simulation > 0 else 0
        report_lines.append(f"  {'Absorbed (Inside)':<30} {self.N_absorbed:<8} {pct_abs_count:<10.2f} {self.E_absorbed:<15.6f} {pct_abs_energy:<10.2f}")
        report_lines.append("-" * 80)
        report_lines.append(f"  {'TOTAL':<30} {self.Nsim:<8} {100.0:<10.2f} {self.E_total_simulation:<15.6f} {100.0:<10.2f}")

        report_lines.append(f"\n{'TRANSMITTED PHOTONS BREAKDOWN':^70}")
        report_lines.append("-" * 80)
        report_lines.append(f"  {'Category':<30} {'Count':<12} {'%':<10} {'Energy (MeV)':<15} {'%':<10}")
        report_lines.append("-" * 80)
        pct_primary_count = (self.N_out_primaries / self.Nsim * 100) if self.Nsim > 0 else 0
        pct_primary_energy = (self.E_out_primaries / self.E_total_simulation * 100) if self.E_total_simulation > 0 else 0
        pct_secondary_count = (self.N_out_secondaries / self.Nsim * 100) if self.Nsim > 0 else 0
        pct_secondary_energy = (self.E_out_secondaries / self.E_total_simulation * 100) if self.E_total_simulation > 0 else 0
        report_lines.append(f"  {'Primaries':<30} {self.N_out_primaries:<12} {pct_primary_count:<10.2f} {self.E_out_primaries:<15.6f} {pct_primary_energy:<10.2f}")
        report_lines.append(f"  {'Secondaries (Build-up)':<30} {self.N_out_secondaries:<12} {pct_secondary_count:<10.2f} {self.E_out_secondaries:<15.6f} {pct_secondary_energy:<10.2f}")
        report_lines.append("-" * 80)
        report_lines.append(f"  {'TOTAL':<30} {self.N_out_total:<12} {pct_out_count:<10.2f} {self.E_out_total:<15.6f} {pct_out_energy:<10.2f}")
              
        # add expected number of non-interacted photons, calculate deviation (in sigma units)
        self.N_out_primaries_theory = self.Nsim * 0.5 ** self.num_hvlX 
        if self.simulation_method == "pdf":
            self.N_out_primaries_theory = self.Nsim * 0.5 ** (self.num_hvlX/self.pef)  # adjust for path extension factor in pdf method
        sigma_primaries_theory = np.sqrt(self.N_out_primaries_theory)
        deviation_primaries = (self.N_out_primaries - self.N_out_primaries_theory) / sigma_primaries_theory if sigma_primaries_theory > 0 else 0
        report_lines.append(f"\n  [Expected primaries:   {int(self.N_out_primaries_theory):<12} ({100*0.5 ** self.num_hvlX:.2f}% of total events)]")
        report_lines.append(f"  [Deviation (sigma units):  {deviation_primaries:<12.2f}]\n")
        # add warning if deviation is larger than 3 sigma (which could indicate a problem with the simulation; 3sigma ~ 99.7% confidence interval for a normal distribution)
        if abs(deviation_primaries) > 3:
            report_lines.append("!  [WARNING: Deviation exceeds 3 sigma! Simulation results may be unreliable.]")

        report_lines.append(f"\n{'ENERGY BALANCE':^70}")
        report_lines.append("-" * 80)
        report_lines.append(f"  Total Simulated Energy:     {self.E_total_simulation:.6f} MeV")
        report_lines.append(f"  Total Detected Energy:      {self.E_processed_simulation:.6f} MeV")
        report_lines.append(f"  Energy Balance Error:       {(self.E_total_simulation - self.E_processed_simulation):.6f} MeV ({((self.E_total_simulation - self.E_processed_simulation) / self.E_total_simulation * 100) if self.E_total_simulation > 0 else 0:.2f} %)")

        report_lines.append(f"\n{'BUILD-UP FACTOR CALCULATION':^70}")
        report_lines.append("-" * 80)
        report_lines.append(f"  Transmitted Energy (E_out_total):      {self.E_out_total:.6f} MeV")
        report_lines.append(f"  Theoretical Energy (primaries; NBG):   {self.E_out_primaries_theory:.6f} MeV")
        report_lines.append(f"  Build-up Factor (B(E)):                {self.buildup_factor:.6f}")
        report_lines.append("=" * 70)

        report_text = "\n".join(report_lines)

        if print_report:
            print(report_text)

        if save_to_file:
            e0_tag = f"{self.E0:.2f}".replace(".", "p")
            timestamp = dt.now().strftime("%m%d_%H%M%S")
            file_path = f"initial_output/run26_{self.simulation_method}_{e0_tag}MeV_{self.num_hvlX}hvl_{timestamp}.txt"

            #* Write to file
            with open(file_path, "w", encoding="utf-8") as f:
                f.write(report_text + "\n")
            if save_to_file:
                print(f"Report saved to: {file_path}")


