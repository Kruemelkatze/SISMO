import pandas as pd
import matplotlib.pyplot as plt
import glob
import numpy as np
import re

# === Pfad zu deinen CSV-Dateien ===
files = glob.glob(
    r"./ACSOS/ACSOS_SISMO/*.csv"
)

print(f"Found {len(files)} CSV files.")

data_by_food = {}

# === CSV-Dateien laden und Mixing Ratio berechnen ===
for filename in files:
    # food=5, food=10, ... aus Dateiname extrahieren
    match = re.search(r"food=(\d+)", filename)

    if not match:
        print(f"Skipped file, no food value found: {filename}")
        continue

    food = int(match.group(1))
    df = pd.read_csv(filename)

    # Prüfen, ob nötige Spalten vorhanden sind
    required_cols = {"blue", "yellow", "green"}
    if not required_cols.issubset(df.columns):
        print(f"Skipped file, missing columns: {filename}")
        print("Columns found:", df.columns.tolist())
        continue

    # Mixing Ratio laut Paper:
    # M(t) = p_green(t) = green / (blue + yellow + green)
    total = df["blue"] + df["yellow"] + df["green"]
    M = df["green"] / total.replace(0, np.nan)
    M = M.fillna(0)

    if food not in data_by_food:
        data_by_food[food] = []

    data_by_food[food].append(M.values)

    print(f"Loaded: {filename} | food={food} | time steps={len(M)}")

print("Food groups found:", sorted(data_by_food.keys()))

if not data_by_food:
    raise ValueError("No valid data found. Check file names and CSV columns.")


# === Hilfsfunktion: Zeitpunkt bis Schwelle ===
def time_to_threshold(M, threshold):
    """
    Gibt den ersten Zeitpunkt zurück, an dem M(t) >= threshold.
    Falls die Schwelle nie erreicht wird, wird NaN zurückgegeben.
    """
    indices = np.where(M >= threshold)[0]
    if len(indices) == 0:
        return np.nan
    return indices[0]


# === Plot: Mean M(t) mit Standardabweichung ===
plt.figure(figsize=(8, 6))

summary_rows = []

for food, runs in sorted(data_by_food.items()):
    # Falls Runs unterschiedliche Länge haben: auf kürzeste Länge kürzen
    min_len = min(len(r) for r in runs)
    runs = np.array([r[:min_len] for r in runs])

    mean_M = np.mean(runs, axis=0)
    std_M = np.std(runs, axis=0)

    time = np.arange(min_len)

    # Mittelwertkurve
    plt.plot(
        time,
        mean_M,
        linewidth=2,
        label=f"{food} food sources"
    )

    # Streuungsband
    plt.fill_between(
        time,
        np.maximum(mean_M - std_M, 0),
        np.minimum(mean_M + std_M, 1),
        alpha=0.2
    )

    # === t0.5 / t0.9 pro Run berechnen ===
    t05_values = []
    t09_values = []
    final_values = []

    for run in runs:
        t05_values.append(time_to_threshold(run, 0.5))
        t09_values.append(time_to_threshold(run, 0.9))
        final_values.append(run[-1])

    t05_values = np.array(t05_values)
    t09_values = np.array(t09_values)
    final_values = np.array(final_values)

    # Zusammenfassung pro Food-Level
    summary_rows.append({
        "Food sources": food,
        "Runs": len(runs),
        "Mean t0.5": np.nanmean(t05_values),
        "Std t0.5": np.nanstd(t05_values),
        "Reached t0.5": np.sum(~np.isnan(t05_values)),
        "Mean t0.9": np.nanmean(t09_values),
        "Std t0.9": np.nanstd(t09_values),
        "Reached t0.9": np.sum(~np.isnan(t09_values)),
        "Mean final M(T)": np.nanmean(final_values),
        "Std final M(T)": np.nanstd(final_values)
    })

plt.xlabel("Time step")
plt.ylabel("Mixing ratio M(t)")
plt.title("Information mixing for different numbers of food sources")
plt.ylim(0, 1)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig("mixing_food_comparison.png", dpi=300)
plt.show()

print("Saved plot as mixing_food_comparison.png")


# === Plot: Small Multiples (Faceted Plot) ===
use_common_x_range = True  # Flag: Use max length across all simulations for x-axis
override_x_range = 50      # If > 0, overrides use_common_x_range and forces this x max

effective_sharex = use_common_x_range or (override_x_range > 0)

num_plots = len(data_by_food)
if num_plots <= 5:
    rows, cols = 1, num_plots
else:
    cols = 3
    rows = (num_plots + cols - 1) // cols

fig, axes = plt.subplots(rows, cols, figsize=(3.5 * cols, 3.5 * rows), sharey=True, sharex=effective_sharex, squeeze=False)
axes = axes.flatten()

global_max_len = 0
if override_x_range > 0:
    global_max_len = override_x_range
elif use_common_x_range:
    global_max_len = max(max(len(r) for r in runs) for runs in data_by_food.values()) - 1

for i, (food, runs) in enumerate(sorted(data_by_food.items())):
    ax = axes[i]
    min_len = min(len(r) for r in runs)
    runs_arr = np.array([r[:min_len] for r in runs])
    
    mean_M = np.mean(runs_arr, axis=0)
    std_M = np.std(runs_arr, axis=0)
    time = np.arange(min_len)
    
    color = f"C{i % 10}"
    ax.plot(time, mean_M, linewidth=2, color=color)
    ax.fill_between(time, np.maximum(mean_M - std_M, 0), np.minimum(mean_M + std_M, 1), alpha=0.2, color=color)
    
    ax.set_title(f"{food} food sources")
    ax.set_ylim(0, 1)
    if override_x_range > 0 or use_common_x_range:
        ax.set_xlim(0, global_max_len)
        
    ax.grid(True, linestyle="--", alpha=0.5)
    
    # X-Achsen-Label für die untersten Plots jeder Spalte
    if i + cols >= num_plots:
        ax.set_xlabel("Time step")
        if override_x_range > 0 or use_common_x_range:
            ax.tick_params(labelbottom=True)
        
    if i % cols == 0:
        ax.set_ylabel("Mixing ratio M(t)")

# Verbleibende leere Subplots verstecken
for j in range(num_plots, len(axes)):
    fig.delaxes(axes[j])

plt.suptitle("Information mixing per food source", y=1.02)
plt.tight_layout()

plt.savefig("mixing_food_comparison_faceted.png", dpi=300, bbox_inches="tight")
plt.show()

print("Saved plot as mixing_food_comparison_faceted.png")


# === Tabelle erstellen ===
summary = pd.DataFrame(summary_rows)

# Schöner runden
summary_rounded = summary.copy()
numeric_cols = summary_rounded.select_dtypes(include=[np.number]).columns
summary_rounded[numeric_cols] = summary_rounded[numeric_cols].round(3)

print("\n=== Summary Table ===")
print(summary_rounded.to_string(index=False))

# Tabelle als CSV speichern
summary_rounded.to_csv("mixing_summary_table.csv", index=False)

print("\nSaved table as mixing_summary_table.csv")


# === Optional: zusätzlicher Plot t0.5 gegen Food Sources ===
plt.figure(figsize=(7, 5))

valid_t05 = summary.dropna(subset=["Mean t0.5"])

plt.errorbar(
    valid_t05["Food sources"],
    valid_t05["Mean t0.5"],
    yerr=valid_t05["Std t0.5"],
    marker="o",
    linewidth=2,
    capsize=4
)

plt.xlabel("Number of food sources")
plt.ylabel("Time to reach M(t) ≥ 0.5")
plt.title("Time to 50% information mixing")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig("time_to_50_mixing.png", dpi=300)
plt.show()

print("Saved plot as time_to_50_mixing.png")

# === Plot: Final Mixing ===
plt.figure(figsize=(7, 5))

plt.errorbar(
    summary["Food sources"],
    summary["Mean final M(T)"],
    yerr=summary["Std final M(T)"],
    marker="o",
    linewidth=2,
    capsize=4
)

plt.xlabel("Number of food sources")
plt.ylabel("Final mixing ratio M(T)")
plt.title("Final information mixing")
plt.ylim(0, 1)
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig("final_mixing_by_food.png", dpi=300)
plt.show()

print("Saved plot as final_mixing_by_food.png")