import pandas as pd
import matplotlib.pyplot as plt
import glob
import numpy as np
import re

# === Pfade ===
sismo_files = glob.glob(
    r"./output/*.color_counts.csv"
)

virus_files = glob.glob(
    r"./output_virus/*.color_counts.csv"
)

print(f"Found {len(sismo_files)} SISMO files.")
print(f"Found {len(virus_files)} Virus files.")


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


# === Mittelwert und Standardabweichung berechnen ===
def mean_std_runs(runs, normalize_time=True):
    if len(runs) == 0:
        raise ValueError("No runs available.")

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
            
    runs = np.array(padded_runs)

    mean_M = np.mean(runs, axis=0)
    std_M = np.std(runs, axis=0)
    
    if normalize_time and max_len > 1:
        time = np.linspace(0, 1, max_len)
    else:
        time = np.arange(max_len)

    return time, mean_M, std_M


# === Daten laden ===
sismo_by_food = load_sismo_by_food(sismo_files)
virus_runs = load_virus_runs(virus_files)

print("SISMO food groups:", sorted(sismo_by_food.keys()))
print("Virus runs:", len(virus_runs))

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

plt.xlabel("Normalized simulation time [0, 1]")
plt.ylabel("Mixing ratio M(t)")
plt.title("Modified Virus Network: information mixing")
plt.ylim(0, 1)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig("virus_mixing.png", dpi=300)
plt.show()

print("Saved virus_mixing.png")


# ============================================================
# 2) SISMO VS VIRUS PLOT
# ============================================================

# Vergleich mit einem SISMO-Food-Level
# Empfehlung: 25, weil dort SISMO zuverlässig mischt
comparison_food = 25

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
    label="Modified Virus Network"
)

plt.fill_between(
    virus_time,
    np.maximum(virus_mean - virus_std, 0),
    np.minimum(virus_mean + virus_std, 1),
    alpha=0.2
)

plt.xlabel("Normalized simulation time [0, 1]")
plt.ylabel("Mixing ratio M(t)")
plt.title("SISMO vs Modified Virus Network")
plt.ylim(0, 1)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig("sismo_vs_virus_mixing.png", dpi=300)
plt.show()

print("Saved sismo_vs_virus_mixing.png")