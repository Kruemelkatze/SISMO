import pandas as pd
import matplotlib.pyplot as plt
import glob
import numpy as np
import re
import os

# === Konfiguration ===
NORMALIZE_TIME = True
TRUNCATE_SERIES = False

# === Pfade ===
virus_file_filter = "nodes=250"
file_id = virus_file_filter + " - "
plot_folder = "plots"

if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)

sismo_files = glob.glob(
    r"./ACSOS/ACSOS_SISMO/*.color_counts.csv"
)

virus_files = glob.glob(
    r"./ACSOS/ACSOS_VIRUS/*" + virus_file_filter + "*.color_counts.csv"
)

print(f"Found {len(sismo_files)} SISMO files.")
print(f"Found {len(virus_files)} Virus files.")

# Vergleich mit einem SISMO-Food-Level
# Empfehlung: 25, weil dort SISMO zuverlässig mischt
comparison_food = 25


# === Hilfsfunktion: M(t) berechnen ===
def compute_mixing(df):
    required_cols = {"blue", "yellow", "green"}

    if not required_cols.issubset(df.columns):
        raise ValueError(f"Missing columns. Found columns: {df.columns.tolist()}")

    total = df["blue"] + df["yellow"] + df["green"]

    # Laut Paper:
    # M(t) = p_green(t) = green / (blue + yellow + green)
    M = df["green"] / total.replace(0, np.nan)
    M = M.fillna(0)

    return M.values


# === SISMO nach Food Sources gruppieren ===
def load_sismo_by_food(files):
    data_by_food = {}

    for filename in files:
        match = re.search(r"food=(\d+)", filename)

        if not match:
            print(f"Skipped SISMO file, no food value found: {filename}")
            continue

        food = int(match.group(1))
        df = pd.read_csv(filename)

        M = compute_mixing(df)

        if food not in data_by_food:
            data_by_food[food] = []

        data_by_food[food].append(M)

    return data_by_food


# === Virus-Runs laden, ohne Food-Gruppierung ===
def load_virus_runs(files):
    runs = []

    for filename in files:
        df = pd.read_csv(filename)
        M = compute_mixing(df)
        runs.append(M)

    return runs


# === Padding/Truncation and Mean/Std Calculation ===
def pad_runs(runs, normalize_time=None, truncate_series=None):
    if normalize_time is None:
        normalize_time = NORMALIZE_TIME
    if truncate_series is None:
        truncate_series = TRUNCATE_SERIES

    if len(runs) == 0:
        raise ValueError("No runs available.")

    if truncate_series:
        # Falls Runs unterschiedliche Länge haben: auf kürzeste Länge kürzen
        min_len = min(len(r) for r in runs)
        runs_arr = np.array([r[:min_len] for r in runs])
        target_len = min_len
    else:
        # Find the maximum length across all runs
        max_len = max(len(r) for r in runs)
        
        # Pad shorter runs with their last value (forward fill)
        padded_runs = []
        for r in runs:
            if len(r) < max_len:
                # Pad with the last value to simulate state persistence
                pad_val = r[-1] if len(r) > 0 else 0
                padded_r = np.pad(r, (0, max_len - len(r)), 'constant', constant_values=pad_val)
                padded_runs.append(padded_r)
            else:
                padded_runs.append(r)
                
        runs_arr = np.array(padded_runs)
        target_len = max_len

    if normalize_time and target_len > 1:
        time = np.linspace(0, 1, target_len)
    else:
        time = np.arange(target_len)

    return time, runs_arr

def mean_std_runs(runs, normalize_time=None, truncate_series=None):
    time, runs_arr = pad_runs(runs, normalize_time, truncate_series)
    mean_M = np.mean(runs_arr, axis=0)
    std_M = np.std(runs_arr, axis=0)
    return time, mean_M, std_M


# === Daten laden ===
sismo_by_food = load_sismo_by_food(sismo_files)
virus_runs = load_virus_runs(virus_files)

print("SISMO food groups:", sorted(sismo_by_food.keys()))
print("Virus runs:", len(virus_runs))
print("SISMO runs:", len(sismo_by_food[comparison_food]))

if not sismo_by_food:
    raise ValueError("No SISMO data found.")

if not virus_runs:
    raise ValueError("No Virus data found.")


# ============================================================
# 1) VIRUS-ONLY PLOT
# ============================================================

virus_time, virus_mean, virus_std = mean_std_runs(virus_runs)

plt.figure(figsize=(8, 6))

plt.plot(
    virus_time,
    virus_mean,
    linewidth=2,
    label="Modified Virus Network"
)

plt.fill_between(
    virus_time,
    np.maximum(virus_mean - virus_std, 0),
    np.minimum(virus_mean + virus_std, 1),
    alpha=0.2
)

plt.xlabel("Normalized simulation time [0, 1]" if NORMALIZE_TIME else "Time step")
plt.ylabel("Mixing ratio M(t)")
plt.title("Modified Virus Network: information mixing")
plt.ylim(0, 1)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig(plot_folder + "/" + file_id + "virus_mixing.png", dpi=300)
plt.show()

print("Saved " + plot_folder + "/" + file_id + "virus_mixing.png")


# ============================================================
# 2) SISMO VS VIRUS PLOT
# ============================================================


if comparison_food not in sismo_by_food:
    raise ValueError(f"SISMO food={comparison_food} not found.")

sismo_time, sismo_mean, sismo_std = mean_std_runs(sismo_by_food[comparison_food])

plt.figure(figsize=(8, 6))

# SISMO
plt.plot(
    sismo_time,
    sismo_mean,
    linewidth=2,
    label=f"SISMO ({comparison_food} food sources)"
)

plt.fill_between(
    sismo_time,
    np.maximum(sismo_mean - sismo_std, 0),
    np.minimum(sismo_mean + sismo_std, 1),
    alpha=0.2
)

# Virus
plt.plot(
    virus_time,
    virus_mean,
    linewidth=2,
    label="Modified Virus Network (250 nodes, 10 sources)"
)

plt.fill_between(
    virus_time,
    np.maximum(virus_mean - virus_std, 0),
    np.minimum(virus_mean + virus_std, 1),
    alpha=0.2
)

plt.xlabel("Normalized simulation time t ∈ [0, 1]" if NORMALIZE_TIME else "Time step")
plt.ylabel("Mixing ratio M(t)")
plt.title("SISMO vs Modified Virus Network")
plt.ylim(0, 1)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig(plot_folder + "/" + file_id + "sismo_vs_virus_mixing.png", dpi=300)
plt.show()

print("Saved " + plot_folder + "/" + file_id + "sismo_vs_virus_mixing.png")


# ============================================================
# 3) VIRUS SPAGHETTI PLOT (Density-based visualization)
# ============================================================

virus_time, virus_runs_arr = pad_runs(virus_runs)

plt.figure(figsize=(8, 6))

# Plot all individual runs with low opacity (spaghetti plot)
for r in virus_runs_arr:
    plt.plot(virus_time, r, color='steelblue', alpha=0.15, linewidth=1.5)

# Calculate and plot the median to show the central tendency
# (Median is less affected by outliers than mean)
virus_median = np.median(virus_runs_arr, axis=0)
plt.plot(virus_time, virus_median, color='navy', linewidth=2, label="Median")

plt.xlabel("Normalized simulation time [0, 1]" if NORMALIZE_TIME else "Time step")
plt.ylabel("Mixing ratio M(t)")
plt.title("Modified Virus Network: Spaghetti Plot (Density)")
plt.ylim(0, 1)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig(plot_folder + "/" + file_id + "virus_mixing_spaghetti.png", dpi=300)
plt.show()

print("Saved " + plot_folder + "/" + file_id + "virus_mixing_spaghetti.png")


# ============================================================
# 4) SISMO VS VIRUS SPAGHETTI PLOT (Density-based visualization)
# ============================================================

sismo_time_arr, sismo_runs_arr = pad_runs(sismo_by_food[comparison_food])

plt.figure(figsize=(8, 6))

# Plot all individual SISMO runs with low opacity
for r in sismo_runs_arr:
    plt.plot(sismo_time_arr, r, color='darkorange', alpha=0.15, linewidth=1.5)

# Plot all individual Virus runs with low opacity
for r in virus_runs_arr:
    plt.plot(virus_time, r, color='steelblue', alpha=0.15, linewidth=1.5)

# Calculate and plot the medians
sismo_median = np.median(sismo_runs_arr, axis=0)
plt.plot(sismo_time_arr, sismo_median, color='darkred', linewidth=2, label=f"SISMO Median ({comparison_food} foods)")
plt.plot(virus_time, virus_median, color='navy', linewidth=2, label="Virus Median")

plt.xlabel("Normalized simulation time [0, 1]" if NORMALIZE_TIME else "Time step")
plt.ylabel("Mixing ratio M(t)")
plt.title("SISMO vs Modified Virus Network: Spaghetti Plot")
plt.ylim(0, 1)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig(plot_folder + "/" + file_id + "sismo_vs_virus_mixing_spaghetti.png", dpi=300)
plt.show()

print("Saved " + plot_folder + "/" + file_id + "sismo_vs_virus_mixing_spaghetti.png")