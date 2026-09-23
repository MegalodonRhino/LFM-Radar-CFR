/*
 * LFM Pulse Compression Radar Simulation with CA-CFAR Detection
 * ---------------------------------------------------------------
 * Simulates a linear-frequency-modulated (chirp) radar pulse, places
 * synthetic targets in range, adds noise, applies a matched filter
 * (pulse compression), then runs Cell-Averaging CFAR to detect targets.
 *
 * Output: CSV file with per-range-bin matched-filter magnitude,
 * CFAR threshold, and detection flag, for the Python plotting script.
 *
 * Build:  gcc -O2 -o lfm_radar lfm_radar.c -lm
 * Run:    ./lfm_radar > radar_output.csv
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define PI 3.14159265358979323846
#define C_LIGHT 299792458.0

typedef struct {
    double fs;            /* sampling frequency (Hz) */
    double pulse_width;   /* chirp duration T (s) */
    double bandwidth;     /* chirp bandwidth B (Hz) */
    int n_samples;        /* total length of simulated receive window */
    double noise_power;   /* noise variance (total, both I & Q) */
} RadarParams;

typedef struct {
    double range_m;
    double amplitude;
} Target;

/* ---------- Gaussian random number via Box-Muller ---------- */
static double gauss_rand(void) {
    double u1 = (rand() + 1.0) / (RAND_MAX + 2.0);
    double u2 = (rand() + 1.0) / (RAND_MAX + 2.0);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
}

/* ---------- Generate LFM chirp: s[n] = exp(j*pi*K*t^2) ---------- */
static void gen_chirp(double complex *chirp, int Np, double fs, double bw, double pw) {
    double K = bw / pw;
    for (int n = 0; n < Np; n++) {
        double t = (double)n / fs - pw / 2.0;
        chirp[n] = cexp(I * PI * K * t * t);
    }
}

/* ---------- Build received signal: delayed target echoes + noise ---------- */
static void gen_received(double complex *rx, int N, const double complex *chirp, int Np,
                          const Target *targets, int n_targets, const RadarParams *p) {
    for (int n = 0; n < N; n++) rx[n] = 0.0;

    for (int t = 0; t < n_targets; t++) {
        int delay = (int)round(2.0 * targets[t].range_m / C_LIGHT * p->fs);
        for (int n = 0; n < Np; n++) {
            int idx = delay + n;
            if (idx >= 0 && idx < N)
                rx[idx] += targets[t].amplitude * chirp[n];
        }
    }

    double sigma = sqrt(p->noise_power / 2.0); /* split across I and Q */
    for (int n = 0; n < N; n++)
        rx[n] += sigma * gauss_rand() + I * sigma * gauss_rand();
}

/*
 * Matched filter via direct convolution with the time-reversed
 * conjugate of the chirp: h[n] = conj(chirp[Np-1-n]).
 * y[k] = sum_m rx[m] * h[k-m],  output length = N + Np - 1
 */
static void matched_filter(const double complex *rx, int N, const double complex *chirp, int Np,
                            double *out_mag, int out_len) {
    for (int k = 0; k < out_len; k++) {
        double complex acc = 0.0;
        int m_lo = (k - Np + 1) > 0 ? (k - Np + 1) : 0;
        int m_hi = (k < N - 1) ? k : (N - 1);
        for (int m = m_lo; m <= m_hi; m++) {
            int hn = k - m;                 /* index into h, 0..Np-1 */
            double complex h = conj(chirp[Np - 1 - hn]);
            acc += rx[m] * h;
        }
        out_mag[k] = cabs(acc);
    }
}

/*
 * CA-CFAR on the power (magnitude^2) profile.
 * guard: guard cells on each side of the cell under test (CUT)
 * train: training/reference cells on each side, beyond the guard cells
 * pfa:   desired probability of false alarm
 */
static void ca_cfar(const double *mag, int N, int guard, int train, double pfa,
                     double *threshold, int *detected) {
    int n_ref = 2 * train;
    double alpha = n_ref * (pow(pfa, -1.0 / n_ref) - 1.0); /* CFAR threshold factor */

    for (int k = 0; k < N; k++) {
        int lead_start = k - guard - train;
        int lead_end   = k - guard - 1;
        int lag_start  = k + guard + 1;
        int lag_end    = k + guard + train;

        double sum = 0.0;
        int count = 0;

        for (int i = lead_start; i <= lead_end; i++) {
            if (i >= 0 && i < N) { sum += mag[i] * mag[i]; count++; }
        }
        for (int i = lag_start; i <= lag_end; i++) {
            if (i >= 0 && i < N) { sum += mag[i] * mag[i]; count++; }
        }

        if (count == 0) {
            threshold[k] = 0.0;
            detected[k] = 0;
            continue;
        }

        double noise_est = sum / count;
        double thresh_power = alpha * noise_est;
        threshold[k] = sqrt(thresh_power); /* back to magnitude domain for plotting */
        detected[k] = (mag[k] * mag[k] > thresh_power) ? 1 : 0;
    }
}

int main(void) {
    srand((unsigned int)time(NULL));

    RadarParams p = {
        .fs = 10e6,            /* 10 MHz sample rate */
        .pulse_width = 10e-6,  /* 10 us chirp */
        .bandwidth = 5e6,      /* 5 MHz bandwidth -> range resolution c/(2B) ~= 30 m */
        .n_samples = 2048,     /* receive window length */
        .noise_power = 0.6     /* noise variance */
    };

    int Np = (int)round(p.pulse_width * p.fs); /* samples per chirp */
    double complex *chirp = malloc(sizeof(double complex) * Np);
    gen_chirp(chirp, Np, p.fs, p.bandwidth, p.pulse_width);

    /* Define synthetic targets: range (m), amplitude */
    Target targets[] = {
        {1500.0, 1.0},
        {1560.0, 0.6},   /* close to target 1: tests range resolution */
        {4200.0, 1.4},
        {7800.0, 0.35}   /* weak target: tests CFAR sensitivity */
    };
    int n_targets = sizeof(targets) / sizeof(targets[0]);

    double complex *rx = malloc(sizeof(double complex) * p.n_samples);
    gen_received(rx, p.n_samples, chirp, Np, targets, n_targets, &p);

    int out_len = p.n_samples + Np - 1;
    double *mag = malloc(sizeof(double) * out_len);
    matched_filter(rx, p.n_samples, chirp, Np, mag, out_len);

    double *threshold = malloc(sizeof(double) * out_len);
    int *detected = malloc(sizeof(int) * out_len);
    ca_cfar(mag, out_len, /*guard=*/4, /*train=*/16, /*pfa=*/1e-4, threshold, detected);

    /*
     * The "full" convolution used by matched_filter() introduces a group
     * delay of (Np-1) samples: a target echo starting at sample `delay`
     * produces its compressed peak at output index `delay + Np - 1`.
     * Shift the index back by (Np-1) so bin/range line up with the
     * original receive window (and drop the leading ramp-up transient
     * where the filter hasn't fully overlapped the echo yet).
     */
    printf("bin,range_m,mf_magnitude,cfar_threshold,detected\n");
    for (int k = Np - 1; k < out_len; k++) {
        int bin = k - (Np - 1);
        double range_m = (double)bin / p.fs * C_LIGHT / 2.0;
        printf("%d,%.4f,%.6f,%.6f,%d\n", bin, range_m, mag[k], threshold[k], detected[k]);
    }

    free(chirp); free(rx); free(mag); free(threshold); free(detected);
    return 0;
}
