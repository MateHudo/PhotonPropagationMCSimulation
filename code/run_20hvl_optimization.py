#import numpy as np
import os
from pathlib import Path
import pandas as pd
from datetime import datetime as dt

def save_row(row_dict, filepath):
    df = pd.DataFrame([row_dict])
    if not os.path.exists(filepath):
        df.to_csv(filepath, index=False)
    else:
        df.to_csv(filepath, mode="a", header=False, index=False)

#* filename and filepath; "where to save results"
filename = f"res_20hvl_{dt.now().strftime('%Y%m%d_%H%M%S')}.csv"
relative_path_to_file = [
    "initial_output/simulation_results", # old path
    "../data/raw",
    "data", #note: use this one, rename and then move it to /raw folder
] [2]
savefile_path = relative_path_to_file + f"/{filename}"
# Stop if file already exists (omitting overwriting/appending to existing file from previous simulations); not relevant 
# if we use unique filename with timestamp.
if Path(savefile_path).exists():
    print(f"WARNING: File already exists: {savefile_path}")
    print("Choose a new filename to avoid overwriting/appending.")
    raise FileExistsError(f"File already exists: {savefile_path}")


#import importlib
from code.src import PhotonBoxPropagationSimulator as Simulator
from code.src import LinearAttenuationCoeff_load as LAC_loader

# load LAC data
fp_mac = "mac_NistXcom/mac_lead.txt"
lead_density = 11.34 # g/cm^3
lac_loader = LAC_loader.LACLoader(fp_mac,lead_density,"Svinec")


#$ simulation parameters
#* I) main
Nsim = 10000000
N_reps = 10

# all energies!
E0_list = [0.1,0.25,0.5,1.0,2.0,5.0,10.0]
#E0_list = [1]

# 20 HVL
N_hvl_list = [20]
n_hvl_y = 60
n_hvl_z = 60

#? all methods? -> No!
#* - pdf and combine are best!
#* - normal and forcing are worse, so ommit them
methods = ["combine", "pdf_man"]

#$ simulate for different pef and weight_min
#! I already have results for pef=2 and w_min = 0.01, so no need for this combination -> but hard to ommit?
path_extension_factor_list = [2,3,4,5]
weight_min_list = [0.01,0.001,0.0001]
#? also for survival_probability?
#survival_probability_list = 0.1

#* II) others
verbose = False
force_first_interaction = True
survival_probability = 0.1



#################
#####  RUN  #####
#################
start_time = dt.now()
total_simulations = len(methods)*len(E0_list)*len(N_hvl_list)*N_reps*len(path_extension_factor_list)*len(weight_min_list) 
sim_idx = - 1
print("\n------------- STARTING SIMULATIONS -------------")
print(f"Total simulations to run: {total_simulations}")
print(f"Simulation parameters:")
print(f"\tMethods: {methods}")
print(f"\tE0 values: {E0_list}")
print(f"\tN_hvl values: {N_hvl_list}")
print(f"\tRepetitions: {N_reps}")
print(f"\tNsim per simulation: {Nsim}")
print("----------------------------------------------\n")
for method in methods:
    for E0 in E0_list:
        for N_hvl in N_hvl_list:
            for rep in range(N_reps):
                for path_extension_factor in path_extension_factor_list:
                    for weight_min in weight_min_list:

                        sim_idx += 1
                        # simulate for choosen setup and gather results
                        print(f"\n\n")
                        print("="*75)
                        print(f"\t\tMethod: {method}, E0: {E0} MeV, N_hvl: {N_hvl}, Rep: {rep+1}/{N_reps}, PEF: {path_extension_factor}, w_min: {weight_min}")
                        time = dt.now() - start_time
                        print(f"\t\t(N_sim_completed={sim_idx}/{total_simulations}; total_time={time.total_seconds()/60:.2f} min)")
                        print("="*75)

                        # ---------------------------------
                        # initialize simulator class instance
                        # ---------------------------------
                        simulator = Simulator.photon_box_propagation_simulator((N_hvl, n_hvl_y, n_hvl_z), lac_loader, E0)
                        config = {
                            "simulation_method": method,
                            "Nsim": Nsim,
                            "Srep_index": rep,
                            "verbose": verbose,
                            "progress_report_frequency": 10,
                            "path_extension_factor": path_extension_factor,
                            "force_first_interaction": force_first_interaction,
                            "weight_min": weight_min,
                            "start_weight": 1.0,
                            "survival_probability": survival_probability
                        }

                        result = simulator.run(config)
                        
                        #! Saving single simulation results
                        #* I) normal storing and saving later --> if used, later needs saving: df=pd.DataFrame(results) and df.to_csv(...)
                        #results.append(result)
                        #* II) real time partial results saving --> safer as intermediate results are saved, not lost if simulation breaks
                        save_row(result, savefile_path)



print("\n\n\n------------- SIMULATIONS COMPLETED -------------")
#note: tole se da izračunat tut v data_analysis!
total_time = dt.now() - start_time
print(f"Total simulation time: {total_time.total_seconds()/60:.2f} min")