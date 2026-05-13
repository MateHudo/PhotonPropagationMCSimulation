# import
import numpy as np
import os
import pandas as pd
from datetime import datetime as dt

# ----------------------------
#  DEFINE INPUT FILEPATH AND NAME
# -----------------------------    
# simulation raw data to process and save
input_filepath = "../data/raw"
input_filename = "remote/res_2to3hvl.csv"
# --------------------
# DEFINE OUTPUT FILEPATH AND NAME
# --------------------------
output_filepath = "../data/processed/"
#output_filename = f"processed_{dt.now().strftime('%Y%m%d_%H%M%S')}.csv"
output_filename = f"processed_2to3hvl.csv"
#$ I would like to check if the output file already exists, and if so: raise an error to avoid overwriting existing data
if os.path.exists(os.path.join(output_filepath, output_filename)):
    raise FileExistsError(f"Output file {output_filename} already exists in {output_filepath}. Please change the output filename to avoid overwriting existing data.")


# ----------------------------
#  LOAD SIMULATION RAW DATA
# ----------------------------
# dataframe with raw data
print(f"\n\nLoading raw simulation data from {os.path.join(input_filepath, input_filename)}")
df = pd.read_csv(f"{input_filepath}/{input_filename}")


# -------------------
#  PROCESS RAW DATA
# -------------------
#  
# cleanup - remove unnecessary columns from the dataframe - based on keys names
columns_to_drop = [
    "box_shape_hvl", "Nsim", "N_steps", "force_first_interaction", "path_extension_factor", "weight_min",
    "survival_probability", "box_size_cm", "verbosity", "Emin_terminate", "N_absorbed", "E_absorbed",
    "N_stps_per_sec", "n_hvl_y", "n_hvl_z", "N_out_total", "E_out_total", "E_out_primaries",
    "E_out_secondaries", "N_backscattered", "E_backscattered", "N_leakage", "E_leakage",

    #"N_out_primaries", #"N_out_secondaries", #"N_stps_per_Nsim", #"Nsim_per_sec"
]
df = df.drop(columns=columns_to_drop, errors="ignore")

# Main aggregated table
summary = (
    df.groupby(["method", "E0", "n_hvl_x"], as_index=False)
      .agg(
          n_rep=("Srep_index", "nunique"),
          t_mean=("simulation_time_sec", "mean"),
          t_std=("simulation_time_sec", "std"),
          t_min=("simulation_time_sec", "min"),
          t_max=("simulation_time_sec", "max"),
          # flag if relative time spread > 5%
          t_unstable=("simulation_time_sec", lambda s: (s.max() - s.min()) / s.mean() > 0.05),
          B_mean=("buildup_factor", "mean"),
          B_std=("buildup_factor", "std"),
          # rounded mean number of escaping secondaries (integer)
          Nsec_mean=("N_out_secondaries", lambda s: int(np.rint(s.mean()))),
          Nsec_std=("N_out_secondaries", "std"),
          Nprim_mean=("N_out_primaries", "mean"),
          Nprim_std=("N_out_primaries", "std"),
          # steps per simulated photon
          stps_per_ph=("N_stps_per_Nsim", "mean"),
          # simulated photons per second
          Nsim_per_sec=("Nsim_per_sec", "mean"),
      )
)

# change "buildup" method name to "normal" and combined_altered/combined to "combine"
summary["method"] = summary["method"].replace({"buildup": "normal"})
summary["method"] = summary["method"].replace({"combined_altered": "combine"}) # combined_altered
summary["method"] = summary["method"].replace({"combined": "combine"}) # combined
# change column n_hlv_x to n_hvl
summary = summary.rename(columns={"n_hvl_x": "n_hvl"})

# calculate FOM = 1 / (sigma_B^2 * t_mean)
#summary["FOM"] = 1.0 / (summary["B_std"]**2 * summary["t_mean"]) # using t_mean
summary["FOM"] = 1.0 / (summary["B_std"]**2 * summary["t_min"]) # using t_min - better! (as it is more robust to outliers in time measurements)

# sort by E0, n_hvl, method
summary = summary.sort_values(["E0", "n_hvl", "method"]).reset_index(drop=True)

# further columns dropping:  n_rep, t_mean, t_std, t_max, t_unstable
summary = summary.drop(columns=[
    "n_rep",
    "t_mean",
    "t_std",
    "t_max",
    "t_unstable",
])
# rename: t_min --> t_sim
summary = summary.rename(columns={"t_min": "t_sim"})


print("\nData processing complete.")

# -- -------------------
# SAVE PROCESSED SUMMARY TABLE TO CSV
# -------------------
summary.to_csv(os.path.join(output_filepath, output_filename), index=False)
print(f"\nProcessed summary table saved to {os.path.join(output_filepath, output_filename)}.")