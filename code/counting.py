
 #! ...
self.N_steps = 0  # photon step is defined as one photon discrete move (or attempted move)
self.N_steps_successful = 0  # N_steps - N_terminated

# out
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

# box - absorbed by the medium or termination based on energy or weight threshold
self.N_box = 0  # number of photons absorbed in the medium (photoeffect or termination, based on w_min or E_min threshold)
self.N_absorbed = 0  # number of photons absorbed in the medium (photoeffect only)
self.E_absorbed = 0.0
self.N_weight_terminated = 0  # number of photons terminated by Russian roulette
self.N_energy_terminated = 0  # number of photons terminated by energy threshold
#? those two below??
self.E_weight_terminated = 0.0  # energy of photons terminated by Russian roulette
self.E_energy_terminated = 0.0  # energy of photons terminated by energy threshold
# counting interactions
self.N_interactions_total = 0
self.N_compton = 0
self.N_photoeffect = 0
self.N_pair_production = 0