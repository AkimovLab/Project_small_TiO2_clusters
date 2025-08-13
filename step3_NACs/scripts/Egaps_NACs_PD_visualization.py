import numpy as np
import matplotlib.pyplot as plt


# === Give full paths to your combined .txt files ===
FILES = {
    100: "../key_outputs/gap_01_100K.txt",
    300: "../key_outputs/gap_01_300K.txt"
}

SYSTEM_IDS = [1, 2, 3, 4, 5, 6,7,8]   # four systems to plot
YLABEL = "PD (1/eV)"
XLABEL = "Energy gap (eV)"
LOG_Y = True
XLIM = (0.0, 5.0)           # adjust if needed

# Load each file into memory
combined_data = {}
for T, filepath in FILES.items():
    data = np.loadtxt(filepath, comments="#", dtype=str)
    combined_data[T] = {
        "system_id": data[:, 1].astype(int),
        "system_label": data[:, 2],
        "bin": data[:, 3].astype(float),
        "dens": data[:, 4].astype(float),
        "cum": data[:, 5].astype(float),
    }

# Make figure (4 rows × 2 columns)
fig, axes = plt.subplots(len(SYSTEM_IDS), len(FILES), figsize=(10, 10), sharex=True, sharey=True)

for i, sid in enumerate(SYSTEM_IDS):
    for j, T in enumerate(FILES.keys()):
        ax = axes[i, j]
        d = combined_data[T]
        mask = d["system_id"] == sid
        if not mask.any():
            ax.text(0.5, 0.5, f"No data\nsystem {sid}",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            continue

        x = d["bin"][mask]
        y = d["dens"][mask]
        ax.plot(x, y)

        if LOG_Y:
            ax.set_yscale("log")
        if XLIM:
            ax.set_xlim(*XLIM)

        if j == 0:
            ax.set_ylabel(YLABEL)
        if i == len(SYSTEM_IDS) - 1:
            ax.set_xlabel(XLABEL)

        ax.set_title(f"System {sid} — {T}K", fontsize=10)

fig.tight_layout()
plt.show()




