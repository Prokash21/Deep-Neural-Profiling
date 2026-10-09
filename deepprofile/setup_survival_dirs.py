import os

base = "SURVIVAL_ANALYSIS"
dirs = [
    os.path.join(base, "data"),
    os.path.join(base, "scripts"),
    os.path.join(base, "notebooks"),
    os.path.join(base, "output", "km_plots"),
    os.path.join(base, "output", "cox_tables"),
    os.path.join(base, "output", "forest_plots"),
]

for d in dirs:
    os.makedirs(d, exist_ok=True)

print(f"Created the following folders under {base}/:")
for d in dirs:
    print("  " + d)
