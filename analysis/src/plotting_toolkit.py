import numpy as np
import matplotlib.pyplot as plt

def show_setup(summary_df, E0, n_hvl, method=None):
    q = (summary_df["E0"] == E0) & (summary_df["n_hvl"] == n_hvl)
    if method is not None:
        q &= (summary_df["method"] == method)

    out = summary_df.loc[q, ["method", "E0", "n_hvl", "t_sim", "B_mean", "B_std", "FOM"]]
    return out.sort_values("method")




def plot_B_vs_hvl(df, E0, methods_used=["all"], x_shift=0.18):
    sub = df[df["E0"] == E0].copy()

    if methods_used != ["all"]:
        sub = sub[sub["method"].isin(methods_used)]
        method_order = [m for m in methods_used if m in sub["method"].unique()]
    else:
        method_order = sorted(sub["method"].unique())

    n_methods = len(method_order)
    if n_methods == 0:
        return
    offsets = np.linspace(-x_shift, x_shift, n_methods) if n_methods > 1 else [0.0]

    plt.figure(figsize=(7, 5))
    for m, dx in zip(method_order, offsets):
        g = sub[sub["method"] == m].sort_values("n_hvl")
        x = g["n_hvl"] + dx
        plt.errorbar(
            x, g["B_mean"], yerr=g["B_std"],
            marker="o", capsize=3, linestyle="", label=m
        )

    plt.xlabel("N_hvl")
    plt.ylabel("Buildup faktor (B)")
    plt.title(f"B(N_hvl) [E0={E0} MeV]")
    plt.legend(fontsize="small", loc="best")
    plt.tight_layout()
    plt.show()

def plot_value_vs_E0(summary_df, n_hvl, y_value="t_sim", yscale="log", xscale="linear"):
    sub = summary_df[summary_df["n_hvl"] == n_hvl].copy()

    plt.figure(figsize=(7, 4))
    for m, g in sub.groupby("method"):
        g = g.sort_values("E0")
        plt.plot(g["E0"], g[y_value], marker="o", label=m)

    plt.xlabel("E0")
    plt.ylabel(f"Mean {y_value} [s]")
    plt.title(f"{y_value} vs E0 (n_hvl={n_hvl})")
    plt.legend(fontsize="small", loc="best")
    plt.yscale(yscale)
    plt.xscale(xscale)
    plt.tight_layout()
    plt.show()

def plot_value_vs_hvl(summary_df, E0, y_value = "t_sim",yscale="log"):
    sub = summary_df[summary_df["E0"] == E0].copy()
    plt.figure(figsize=(7, 4))
    for m, g in sub.groupby("method"):
        g = g.sort_values("n_hvl")
        plt.plot(g["n_hvl"], g[y_value], marker="o", label=m)
    plt.xlabel("N_hvl")
    plt.ylabel(f"Mean {y_value} [s]")
    plt.title(f"{y_value} vs N_hvl (E0={E0})")
    plt.legend(fontsize="small", loc="best")
    plt.yscale(yscale)
    plt.tight_layout()
    plt.show()
