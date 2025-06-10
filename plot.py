import json
import os
import matplotlib.pyplot as plt
import sys
import re

# === CONFIGURATION ===
# Directory containing the log files
LOG_DIR = sys.argv[1]
# Metrics to plot
METRICS = ["particles::collide", "FPS"]
LABELS = ["Czas poświęcony na obliczenia kolizji[sek]", "klatki na sekundę[FPS]"]
# File pattern filter (optional)
LOG_FILE_EXTENSION = ".json"

def read_log_file(filepath):
    with open(filepath, "r") as f:
        data = json.load(f)
    return data["measurements"]

def extract_metric_series(measurements, metric):
    return [float(entry[metric]) for entry in measurements if metric in entry]

import re

# Manual swap map: filename suffix -> desired label
SWAP_LABELS = {
    "seq": "omp",
    "omp": "seq"
}

def plot_metrics(log_files, metrics):
    all_measurements = []

    # Load all measurement data
    for log_file in log_files:
        measurements = read_log_file(os.path.join(LOG_DIR, log_file))
        all_measurements.append(measurements)

    # Determine shortest measurement list
    min_len = min(len(m) for m in all_measurements)

    for idx, metric in enumerate(metrics):
        plt.figure(figsize=(10, 5))

        for log_file, measurements in zip(log_files, all_measurements):
            match = re.match(r"(\d+)(.*)\.json", log_file)
            if match:
                num = match.group(1)
                suffix = match.group(2)
                # Swap the label
                label = SWAP_LABELS.get(suffix, suffix)
                title = f"Test {num}"
            else:
                label = log_file.replace(".json", "")
                title = "Test"

            values = extract_metric_series(measurements[:min_len], metric)
            x_vals = [i / 8 for i in range(len(values))]
            plt.plot(x_vals, values, label=label)

        plt.title(title)
        plt.xlabel("Czas symulacji [sek]")
        plt.ylabel(LABELS[idx])
        plt.legend(title="Implementacja")
        plt.grid(True)
        plt.tight_layout()
        plt.show()


def main():
    log_files = [f for f in os.listdir(LOG_DIR) if f.endswith(LOG_FILE_EXTENSION)]
    if not log_files:
        print("No log files found.")
        return
    plot_metrics(log_files, METRICS)

if __name__ == "__main__":
    main()

