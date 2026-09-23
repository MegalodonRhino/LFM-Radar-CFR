# LFM Pulse Compression Radar Simulation with CFAR Target Detection

![Detection Plot](radar_detection_plot.png)

## What this is

I built this project to get hands-on with a core radar signal processing
technique: pulse compression. I write the full pipeline in C — chirp
generation, target/noise simulation, matched filtering, and adaptive
detection — then visualize the results in Python.

The idea I wanted to explore: real radars can't just blast out a short,
high-power pulse to get good range resolution (the hardware can't handle
the peak power). Instead, they transmit a long, low-power pulse whose
frequency sweeps linearly over time — a chirp — and recover the resolution
afterward with a matched filter. I wanted to actually implement that
trade-off myself instead of just reading about it, and then layer on
CFAR detection to see how a receiver decides "target" vs. "noise" without
a human eyeballing a plot.

## How I built it

**1. Chirp generation.** I generate a linear-frequency-modulated (LFM)
pulse, `s(t) = exp(jπKt²)` where `K = B/T`, using my configured pulse
width `T` and bandwidth `B`.

**2. Scene simulation.** I place a handful of synthetic targets at chosen
ranges and amplitudes — including a closely-spaced pair and one weak
target — delay each one's chirp copy by `2R/c` samples, sum them into a
receive window, and add complex Gaussian noise to stand in for thermal
noise.

**3. Matched filtering (pulse compression).** I convolve the received
signal with the time-reversed conjugate of my transmitted chirp. This is
the step that collapses the long, spread-out chirp back into a narrow,
high-SNR peak at each target's true range — the resolution after
compression works out to `c / (2B)`.

**4. CA-CFAR detection.** I implemented Cell-Averaging Constant False
Alarm Rate detection: for every range bin, I average the power in a ring
of surrounding training cells (skipping guard cells right next to the
cell under test so a target's own sidelobes don't bias the estimate),
then set an adaptive threshold `α × noise_estimate`, where `α` comes from
the false-alarm probability I want. This was the part I found most
interesting — a fixed threshold either misses weak targets or drowns in
false alarms once the noise floor shifts, and CFAR sidesteps that by
constantly re-estimating the local noise.

**5. Visualization.** I pipe the C program's CSV output into a Python
script that plots the matched-filter range profile, the CFAR threshold
curve, and the resulting detections, so I can actually see the pipeline
working end to end.

## A bug I ran into (and why I think it's worth mentioning)

My first version had target detections showing up about 1.5 km further
out than where I'd actually placed the targets. It turned out my matched
filter — implemented as a full convolution — has a built-in group delay
of `Np - 1` samples: a target's compressed peak lands at
`delay + (Np - 1)`, not at `delay`. I fixed the range calculation to
account for that shift. It's a classic radar-sim gotcha, and I'm noting
it here because I think working through *why* it happened taught me more
about matched filtering than getting the numbers right on the first try
would have.

## Build & run

```bash
gcc -O2 -Wall -o lfm_radar lfm_radar.c -lm
./lfm_radar > radar_output.csv
python3 visualize.py radar_output.csv
```

Needs `matplotlib` on the Python side (`pip install matplotlib`).

Note: I seed the noise from `time(NULL)`, so detections near the weak
target can flicker slightly between runs — that's expected CFAR behavior
near threshold, not a bug. For a reproducible run I can point to in a
report, I'd swap that for a fixed seed like `srand(42);`.

## Parameters I used

| Parameter | Meaning | Value |
|---|---|---|
| `fs` | Sample rate | 10 MHz |
| `pulse_width` | Chirp duration `T` | 10 µs |
| `bandwidth` | Chirp bandwidth `B` | 5 MHz → range resolution ≈ 30 m |
| `n_samples` | Receive window length | 2048 |
| `noise_power` | Noise variance | 0.6 |
| `targets[]` | Range (m) / amplitude pairs | 4 targets, incl. a close pair and a weak one |
| CFAR `guard`, `train`, `pfa` | Detector tuning | 4, 16, 1e-4 |

## Where I'd take this next

- Swap the direct O(N·Np) convolution for an FFT-based matched filter and
  compare runtime, once I re-derive the delay bookkeeping for that case.
- Sweep `Pfa` across many noise realizations and plot a proper ROC curve
  to characterize the detector.
- Try OS-CFAR or GO-CFAR and show where CA-CFAR breaks down when a second
  target lands inside the training window.
- Add Doppler processing — per-pulse phase shift for moving targets and a
  2D range-Doppler map.
- Apply a window (Taylor/Hamming) to the chirp to suppress range
  sidelobes, and compare before/after.
