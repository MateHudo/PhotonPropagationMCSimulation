 
#$ ###########################################################################
#$ ################ NEW PLOTTING TOOLKIT FOR PLOTTING BUILDUP RESULTS #########
#$ ################### (used in AnalyseSimulationData.ipynb) ##################
#$ ###########################################################################

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

############################################################
# DEFINING METHODS AND STD MAPS, Y LABELS, COLORS, ...
############################################################
METHOD_ORDER = ["normal", "pdf", "forcing", "combine"]
STD_MAP = {
    "B_mean": "B_std",
    "Nsec_mean": "Nsec_std",
    "Nprim_mean": "Nprim_std",
}
# using slovenian names...
Y_LABELS = {
    "B_mean": "Buildup faktor B",
    "FOM": "FOM",
    "t_sim": "Simulacijski čas [s]",
    "Nsec_mean": r"$N_{sek}$",
    "Nprim_mean": r"$N_{prim}$",
    "stps_per_ph": "Koraki na foton",
    "Nsim_per_sec": "Št. simuliranih fotonov na sekundo",
}
# titles - also in slovenian
TITLE_MAP = {
    "B_mean": "Buildup faktor B",
    "FOM": "Figura zasluge (FOM)",
    "t_sim": "Simulacijski čas",
    "Nsec_mean": r"Povprečno število sekundarnih fotonov $N_{sek}$",
    "Nprim_mean": r"Povprečno število primarnih fotonov $N_{prim}$",
    "stps_per_ph": "Koraki na foton",
    "Nsim_per_sec": "Št. simuliranih fotonov na sekundo",
}

# colors #??
METHOD_COLORS = {m: c for m, c in zip(METHOD_ORDER, plt.rcParams["axes.prop_cycle"].by_key()["color"])}


####################################################################################
# SHOW SETUP - method to show the setup for a given E0, n_hvl, and optionally method
################################################################################
def show_setup(df, E0, n_hvl, method=None):
    q = (df["E0"] == E0) & (df["n_hvl"] == n_hvl)
    if method is not None:
        q &= (df["method"] == method)

    out = df.loc[q, ["method", "E0", "n_hvl", "t_sim", "B_mean", "B_std", "FOM"]]
    return out.sort_values("method")


#####################################################################
#### HELP METHODS FOR PLOTTING (used in main plotting functions below)
#####################################################################
def _get_methods(sub_df, methods="all"):
    if methods == "all":
        return [m for m in METHOD_ORDER if m in sub_df["method"].unique()]
    return [m for m in methods if m in sub_df["method"].unique()]


def _savefig(save_path=None, dpi=300):
    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")


def _energy_colors(energies, cmap_name="cividis"):
    cmap = plt.get_cmap(cmap_name)
    vals = np.linspace(0.15, 0.9, len(energies))  # avoid too pale/dark extremes
    return {E: cmap(v) for E, v in zip(sorted(energies), vals)}

def _nhvl_colors(nhvl_values, cmap_name="cividis"):
    cmap = plt.get_cmap(cmap_name)
    vals = np.linspace(0.15, 0.9, len(nhvl_values))
    return {h: cmap(v) for h, v in zip(sorted(nhvl_values), vals)}

def _bw_style(i):
    linestyles = ["-", "--", "-.", ":"]
    markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    return linestyles[i % len(linestyles)], markers[i % len(markers)]

def _centers_to_edges(vals, log=False):
    vals = np.asarray(sorted(vals), dtype=float)
    if vals.size < 2:
        d = 0.5 if not log else vals[0] * 0.2
        return np.array([vals[0] - d, vals[0] + d], dtype=float)

    if log:
        if np.any(vals <= 0):
            raise ValueError("Log edges require positive values.")
        mids = np.sqrt(vals[:-1] * vals[1:])  # geometric means
        left = vals[0]**2 / mids[0]
        right = vals[-1]**2 / mids[-1]
    else:
        mids = 0.5 * (vals[:-1] + vals[1:])
        left = vals[0] - (mids[0] - vals[0])
        right = vals[-1] + (vals[-1] - mids[-1])

    return np.concatenate([[left], mids, [right]])




###########################################################
######## plot metric vs hvl or E0
###########################################################
def plot_metric_vs_hvl(
    df, 
    E0, 
    y="B_mean", 
    methods="all", 
    show_std=True,
    yscale="linear", 
    linestyle="-", 
    marker="o",
    x_offset=True, 
    offset_width=0.18, 
    figsize=(7, 4),
    title=None, 
    save_path=None
):
    sub = df[df["E0"] == E0].copy()
    methods_used = _get_methods(sub, methods)
    if not methods_used:
        raise ValueError("No methods left after filtering.")

    offsets = np.linspace(-offset_width, offset_width, len(methods_used)) if (x_offset and len(methods_used) > 1) else np.zeros(len(methods_used))
    std_col = STD_MAP.get(y) if show_std else None

    fig, ax = plt.subplots(figsize=figsize)
    for m, dx in zip(methods_used, offsets):
        g = sub[sub["method"] == m].sort_values("n_hvl")
        x = g["n_hvl"].to_numpy(dtype=float) + dx
        yv = g[y].to_numpy(dtype=float)
        yerr = g[std_col].to_numpy(dtype=float) if (std_col in g.columns if std_col else False) else None
        ax.errorbar(x, yv, yerr=yerr, marker=marker, linestyle=linestyle, capsize=3, color=METHOD_COLORS[m], label=m)

    ax.set_xlabel("n_hvl")
    ax.set_ylabel(Y_LABELS.get(y, y))
    ax.set_title(title or f"{y} vs n_hvl (E0={E0} MeV)")
    ax.set_yscale(yscale)
    ax.legend(fontsize="small")
    plt.tight_layout()
    _savefig(save_path)
    plt.show()


def plot_metric_vs_E0(
    df, n_hvl, y="B_mean", methods="all", show_std=True,
    yscale="linear", xscale="log", linestyle="-", marker="o",
    figsize=(7, 4), title=None, save_path=None
):
    sub = df[df["n_hvl"] == n_hvl].copy()
    methods_used = _get_methods(sub, methods)
    if not methods_used:
        raise ValueError("No methods left after filtering.")

    std_col = STD_MAP.get(y) if show_std else None

    fig, ax = plt.subplots(figsize=figsize)
    for m in methods_used:
        g = sub[sub["method"] == m].sort_values("E0")
        yerr = g[std_col].to_numpy(dtype=float) if (std_col in g.columns if std_col else False) else None
        ax.errorbar(g["E0"], g[y], yerr=yerr, marker=marker, linestyle=linestyle, capsize=3, color=METHOD_COLORS[m], label=m)

    ax.set_xlabel("E0 [MeV]")
    ax.set_ylabel(Y_LABELS.get(y, y))
    ax.set_title(title or f"{y} vs E0 (n_hvl={n_hvl})")
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.legend(fontsize="small")
    plt.tight_layout()
    _savefig(save_path)
    plt.show()



################################################
## FOM gain vs hvl (relative to baseline method)
################################################
def plot_FOM_gain_vs_hvl(
    df, E0, baseline="normal", methods="all",
    yscale="log", linestyle="-", marker="o",
    x_offset=True, offset_width=0.18, figsize=(7, 4),
    save_path=None
):
    sub = df[df["E0"] == E0].copy()
    base = sub[sub["method"] == baseline][["n_hvl", "FOM"]].rename(columns={"FOM": "FOM_base"})
    sub = sub.merge(base, on="n_hvl", how="left")
    sub["FOM_gain"] = sub["FOM"] / sub["FOM_base"]

    methods_used = _get_methods(sub, methods)
    if not methods_used:
        raise ValueError("No methods left after filtering.")

    offsets = np.linspace(-offset_width, offset_width, len(methods_used)) if (x_offset and len(methods_used) > 1) else np.zeros(len(methods_used))

    fig, ax = plt.subplots(figsize=figsize)
    for m, dx in zip(methods_used, offsets):
        g = sub[sub["method"] == m].sort_values("n_hvl")
        x = g["n_hvl"].to_numpy(dtype=float) + dx
        ax.plot(x, g["FOM_gain"], marker=marker, linestyle=linestyle, color=METHOD_COLORS[m], label=m)

    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    ax.set_xlabel("n_hvl")
    ax.set_ylabel(f"FOM / FOM_{baseline}")
    ax.set_title(f"FOM gain vs n_hvl (E0={E0} MeV, baseline={baseline})")
    ax.set_yscale(yscale)
    ax.legend(fontsize="small")
    plt.tight_layout()
    _savefig(save_path)
    plt.show()




########################################################################
######## MAIN PLOTTING METHODS (called from AnalyseSimulationData.ipynb)
########################################################################

def plot_B_vs_hvl_multi_E0_single_method(
    df,
    method="normal",
    E0_values="all",
    show_std=False,
    yscale="linear",
    linestyle="-",
    marker="o",
    figsize=(7, 4),
    title=None,
    save_path=None,
    cmap_name="cividis",
    grayscale_safe=False,
):
    sub = df[df["method"] == method].copy()
    if E0_values != "all":
        sub = sub[sub["E0"].isin(E0_values)]

    if sub.empty:
        raise ValueError("No data after filtering. Check method and E0_values.")

    energies = sorted(sub["E0"].unique())
    e_colors = _energy_colors(energies, cmap_name=cmap_name)

    fig, ax = plt.subplots(figsize=figsize)

    for i, E0 in enumerate(energies):
        g = sub[sub["E0"] == E0].sort_values("n_hvl")
        yerr = g["B_std"].to_numpy(dtype=float) if (show_std and "B_std" in g.columns) else None

        if grayscale_safe:
            ls_i, mk_i = _bw_style(i)
            color_i = "black"
        else:
            ls_i, mk_i = linestyle, marker
            color_i = e_colors[E0]

        ax.errorbar(
            g["n_hvl"], g["B_mean"], yerr=yerr,
            linestyle=ls_i, marker=mk_i, capsize=3,
            color=color_i, label=f"E0={E0:g} MeV"
        )

    ax.set_xlabel("N_hvl")
    ax.set_ylabel("Buildup faktor B")
    ax.set_title(title or f"B(N_hvl)")
    ax.set_yscale(yscale)
    ax.legend(fontsize="small", ncol=2)
    plt.tight_layout()
    _savefig(save_path)
    plt.show()


def plot_B_vs_E0_multi_nhvl_single_method(
    df,
    method="normal",
    nhvl_values="all",
    show_std=False,
    xscale="log",
    yscale="linear",
    linestyle="-",
    marker="o",
    figsize=(7, 4),
    title=None,
    save_path=None,
    cmap_name="cividis",
    grayscale_safe=False,
):
    """
    Plot B_mean vs E0 for one chosen method.
    Multiple lines correspond to different n_hvl values.
    """
    sub = df[df["method"] == method].copy()
    if nhvl_values != "all":
        sub = sub[sub["n_hvl"].isin(nhvl_values)]

    if sub.empty:
        raise ValueError("No data after filtering. Check method and nhvl_values.")

    hvls = sorted(sub["n_hvl"].unique())
    h_colors = _nhvl_colors(hvls, cmap_name=cmap_name)

    fig, ax = plt.subplots(figsize=figsize)

    for i, h in enumerate(hvls):
        g = sub[sub["n_hvl"] == h].sort_values("E0")
        yerr = g["B_std"].to_numpy(dtype=float) if (show_std and "B_std" in g.columns) else None

        if grayscale_safe:
            ls_i, mk_i = _bw_style(i)
            color_i = "black"
        else:
            ls_i, mk_i = linestyle, marker
            color_i = h_colors[h]

        ax.errorbar(
            g["E0"], g["B_mean"], yerr=yerr,
            linestyle=ls_i, marker=mk_i, capsize=3,
            color=color_i, label=f"n_hvl={h:g}"
        )

    ax.set_xlabel("E0 [MeV]")
    ax.set_ylabel("Buildup factor B")
    ax.set_title(title or f"B vs E0 (method='{method}')")
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.legend(fontsize="small", ncol=2)

    plt.tight_layout()
    _savefig(save_path)
    plt.show()


def plot_B_vs_hvl_single_E0_multiple_methods(
    df,
    E0,
    methods="all",
    show_std=False,
    yscale="linear",
    linestyle="-",
    marker="o",
    figsize=(7, 4),
    title=None,
    save_path=None,
    grayscale_safe=False,
    x_shift=0.18
):
    """
    Plot B_mean vs n_hvl for one chosen energy.
    Multiple lines correspond to different methods.
    """
    sub = df[df["E0"] == E0].copy()

    if sub.empty:
        raise ValueError(f"No data for E0={E0}.")

    methods_used = _get_methods(sub, methods)
    if not methods_used:
        raise ValueError("No methods left after filtering.")

    fig, ax = plt.subplots(figsize=figsize)

    for i, m in enumerate(methods_used):
        g = sub[sub["method"] == m].sort_values("n_hvl")
        yerr = g["B_std"].to_numpy(dtype=float) if (show_std and "B_std" in g.columns) else None

        if grayscale_safe:
            ls_i, mk_i = _bw_style(i)
            color_i = "black"
        else:
            ls_i, mk_i = linestyle, marker
            color_i = METHOD_COLORS[m]

        ax.errorbar(
            g["n_hvl"], g["B_mean"], yerr=yerr,
            linestyle=ls_i, marker=mk_i, capsize=3,
            color=color_i, label=m
        )

    ax.set_xlabel("N_hvl")
    ax.set_ylabel("Buildup factor B")
    ax.set_title(title or f"B vs N_hvl (E0={E0:g} MeV)")
    ax.set_yscale(yscale)
    ax.legend(fontsize="small")
    plt.tight_layout()
    _savefig(save_path)
    plt.show()


def plot_heatmap_E0_nhvl(
    df, 
    method, 
    value="FOM", 
    cmap="viridis",
    log_color=False, 
    annotate=True, 
    figsize=(7, 5),
    save_path=None, 
    log_x=True, 
    show_cell_grid=True
):
    sub = df[df["method"] == method].copy()
    pvt = sub.pivot(index="n_hvl", columns="E0", values=value).sort_index().sort_index(axis=1)

    x_centers = pvt.columns.to_numpy(dtype=float)   # E0
    y_centers = pvt.index.to_numpy(dtype=float)     # n_hvl
    x_edges = _centers_to_edges(x_centers, log=log_x)
    y_edges = _centers_to_edges(y_centers, log=False)

    arr = pvt.to_numpy(dtype=float)
    if log_color:
        arr = np.where(arr > 0, np.log10(arr), np.nan)

    fig, ax = plt.subplots(figsize=figsize)

    # turn off style grid (it usually causes the "off" look on heatmaps)
    ax.grid(False)

    im = ax.pcolormesh(
        x_edges, y_edges, arr,
        cmap=cmap,
        shading="flat",
        edgecolors=("w" if show_cell_grid else "face"),
        linewidth=(0.6 if show_cell_grid else 0.0),
    )

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(f"log10({value})" if log_color else value)

    ax.set_xticks(x_centers)
    ax.set_yticks(y_centers)
    ax.set_xlabel("E0 [MeV]")
    ax.set_ylabel("n_hvl")
    ax.set_title(f"{value}(E0, n_hvl) \nmethod: {method}")

    if log_x:
        ax.set_xscale("log")

    if annotate:
        for i, y in enumerate(y_centers):
            for j, x in enumerate(x_centers):
                v = pvt.iloc[i, j]
                if pd.notna(v):
                    ax.text(x, y, f"{v:.4g}", ha="center", va="center", fontsize=7, color="k")

    plt.tight_layout()
    _savefig(save_path)
    plt.show()


def plot_surface_E0_nhvl(
    df,
    method="combine",
    value="FOM",
    log_x=True,
    log_z=False,
    cmap="viridis",
    elev=28,
    azim=-130,
    figsize=(8, 6),
    save_path=None,
):
    sub = df[df["method"] == method].copy()
    pvt = sub.pivot(index="n_hvl", columns="E0", values=value).sort_index().sort_index(axis=1)

    E = pvt.columns.to_numpy(dtype=float)   # x centers (E0)
    H = pvt.index.to_numpy(dtype=float)     # y centers (n_hvl)
    Z = pvt.to_numpy(dtype=float)

    X, Y = np.meshgrid(E, H)

    X_plot = np.log10(X) if log_x else X
    Z_plot = np.where(Z > 0, np.log10(Z), np.nan) if log_z else Z

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")

    surf = ax.plot_surface(X_plot, Y, Z_plot, cmap=cmap, linewidth=0, antialiased=True, alpha=0.95)

    cbar = fig.colorbar(surf, ax=ax, shrink=0.72, pad=0.08)
    cbar.set_label(f"log10({value})" if log_z else value)

    if log_x:
        xt = np.log10(E)
        ax.set_xticks(xt)
        ax.set_xticklabels([f"{v:g}" for v in E])
        ax.set_xlabel("E0 [MeV] (log axis)")
    else:
        ax.set_xlabel("E0 [MeV]")

    ax.set_ylabel("n_hvl")
    ax.set_zlabel(f"log10({value})" if log_z else value)
    ax.set_title(f"3D surface: {value}(E0, n_hvl), method={method}")

    ax.view_init(elev=elev, azim=azim)
    plt.tight_layout()
    _savefig(save_path)
    plt.show()