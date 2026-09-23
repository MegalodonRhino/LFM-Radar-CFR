"""
Visualize LFM pulse-compression radar output + CA-CFAR detections.

Usage:
    ./lfm_radar > radar_output.csv
    python3 visualize.py radar_output.csv
"""

import sys
import csv
import matplotlib.pyplot as plt


def load_csv(path):
    bins, ranges, mag, thresh, det = [], [], [], [], []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            bins.append(int(row["bin"]))
            ranges.append(float(row["range_m"]))
            mag.append(float(row["mf_magnitude"]))
            thresh.append(float(row["cfar_threshold"]))
            det.append(int(row["detected"]))
    return bins, ranges, mag, thresh, det


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "radar_output.csv"
    bins, ranges, mag, thresh, det = load_csv(path)

    ranges_km = [r / 1000.0 for r in ranges]
    det_ranges = [r for r, d in zip(ranges_km, det) if d]
    det_mags = [m for m, d in zip(mag, det) if d]

    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)

    # --- Top: matched filter output + CFAR threshold + detections ---
    ax = axes[0]
    ax.plot(ranges_km, mag, color="#1f77b4", linewidth=0.8, label="Matched filter output")
    ax.plot(ranges_km, thresh, color="#d62728", linewidth=0.8,
             linestyle="--", label="CFAR adaptive threshold")
    ax.scatter(det_ranges, det_mags, color="#2ca02c", s=25, zorder=5, label="Detections")
    ax.set_ylabel("Magnitude")
    ax.set_title("LFM Pulse Compression Radar — Range Profile & CA-CFAR Detection")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)

    # --- Bottom: detection flag as a stem/impulse plot ---
    ax2 = axes[1]
    ax2.vlines(det_ranges, 0, 1, color="#2ca02c", linewidth=1.5)
    ax2.set_ylim(0, 1.2)
    ax2.set_yticks([0, 1])
    ax2.set_ylabel("Detected")
    ax2.set_xlabel("Range (km)")
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    out_path = "radar_detection_plot.png"
    plt.savefig(out_path, dpi=150)
    print(f"Saved plot to {out_path}")
    print(f"Detections at ranges (km): {[round(r, 3) for r in det_ranges]}")


if __name__ == "__main__":
    main()
