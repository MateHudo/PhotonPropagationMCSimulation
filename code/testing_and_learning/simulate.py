import numpy as np

#import importlib
from src import PhotonBoxPropagationSimulator as Simulator
from src import LinearAttenuationCoeff_load as LAC_loader

# lac
fp_mac = "../mac_NistXcom/mac_lead.txt"
lead_density = 11.34 # g/cm^3
lac_loader = LAC_loader.LACLoader(fp_mac,lead_density,"Svinec")
print(f"\n----------------------- STARTING SIMULATION -----------------------\n")
#print(f"Loaded LAC data for {lac_loader.material_name} with density {lead_density} g/cm^3.")

#* define parameters:
# side lengths - note: photons start in the middle, so they need to travel only half of the box length to reach the edge!
n_hvl_Y = 20
n_hvl_Z = 20
# box dimensions in HVL units
n_hvl_X = 6
# energy of incident photons
E0 = 1 # MeV
simulation_method = ["simulate",  "buildup", "pdf", "forcing", "combined"] [3] #! change the index to select the method    
Nsim = 100000
verbose = False
path_extension_factor = 2
force_first_interaction = True
weight_min = 0.01
survival_probability = 0.1
#? want reproducible results? set the seed for the random number generator
#np.random.seed(1234355)

#* initialize simulator class instance
simulator = Simulator.photon_box_propagation_simulator((n_hvl_X, n_hvl_Y, n_hvl_Z), lac_loader, E0)

config = {
    "simulation_method": simulation_method,
    "Nsim": Nsim,
    "Srep_index": 0,
    "verbose": verbose,
    "progress_report_frequency": 10,
    "path_extension_factor": path_extension_factor,
    "force_first_interaction": force_first_interaction,
    "weight_min": weight_min,
    "start_weight": 1.0,
    "survival_probability": survival_probability
}
print(f"\nconfig parameters:\n")
for key, value in config.items():
    print(f"{key}: {value}")


# --------------------------------
# INITIALIZE SIMULATION AND RUN IT
# --------------------------------
simulator = Simulator.photon_box_propagation_simulator((n_hvl_X, n_hvl_Y, n_hvl_Z), lac_loader, E0)
sim_info = simulator.run(config)

print("")
for key, value in sim_info.items():
    print(f"{key}: {value}")




