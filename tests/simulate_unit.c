#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../engine/fdtd/simulate.h"

int main(void) {
    fdtd_sim_config_t cfg = {
        .nx = 50,
        .ny = 50,
        .nz = 50,
        .dx = 1e-3,
        .courant_safety = 0.99,
        .n_steps = 200,
        .npml = 8,
        .src_i = 25,
        .src_j = 25,
        .src_k = 25,
        .src_component = FIELD_EZ,
        .t0_over_dt = 40.0,
        .tau_over_dt = 8.0,
        .probe_i = 31,
        .probe_j = 25,
        .probe_k = 25,
        .probe_component = FIELD_EZ
    };

    int result = 1;
    double fs_a = 0.0;
    double fs_b = 0.0;
    double peak = 0.0;
    double early_peak = 0.0;
    double min_value = 0.0;
    double max_value = 0.0;

    double *a = calloc(cfg.n_steps, sizeof(*a));
    double *b = calloc(cfg.n_steps, sizeof(*b));

    if (!a || !b) {
        fprintf(stderr, "FAIL: could not allocate test buffers\n");
        goto cleanup;
    }

    if (fdtd_simulate_echo(&cfg, a, &fs_a) != 0) {
        fprintf(stderr, "FAIL: first simulation failed\n");
        goto cleanup;
    }

    if (fdtd_simulate_echo(&cfg, b, &fs_b) != 0) {
        fprintf(stderr, "FAIL: second simulation failed\n");
        goto cleanup;
    }

    if (!isfinite(fs_a) || fs_a <= 0.0 || !isfinite(fs_b) || fs_b <= 0.0) {
        fprintf(stderr, "FAIL: invalid sample rates: %.9e and %.9e\n", fs_a, fs_b);
        goto cleanup;
    }

    if (fs_a != fs_b) {
        fprintf(stderr, "FAIL: sample rate changed between identical runs\n");
        goto cleanup;
    }

    min_value = a[0];
    max_value = a[0];

    for (size_t i = 0; i < cfg.n_steps; i++) {
        if (!isfinite(a[i]) || !isfinite(b[i])) {
            fprintf(stderr, "FAIL: non-finite sample at index %zu\n", i);
            goto cleanup;
        }

        if (a[i] != b[i]) {
            fprintf(stderr, "FAIL: nondeterministic sample at index %zu: %.17g vs %.17g\n", i, a[i], b[i]);
            goto cleanup;
        }

        double magnitude = fabs(a[i]);

        if (magnitude > peak) peak = magnitude;
        if (i < 20 && magnitude > early_peak) early_peak = magnitude;
        if (a[i] < min_value) min_value = a[i];
        if (a[i] > max_value) max_value = a[i];
    }

    if (peak <= 1e-9) {
        fprintf(stderr, "FAIL: probe did not receive a meaningful signal; peak=%.9e\n", peak);
        goto cleanup;
    }

    if (max_value - min_value <= 1e-9) {
        fprintf(stderr, "FAIL: probe output is effectively constant\n");
        goto cleanup;
    }

    if (early_peak >= peak * 1e-3) {
        fprintf(stderr, "FAIL: early signal exceeds this test's tolerance; early=%.9e, peak=%.9e\n", early_peak, peak);
        goto cleanup;
    }

    if (memcmp(a, b, cfg.n_steps * sizeof(*a)) != 0) {
        fprintf(stderr, "FAIL: identical runs are not bit-for-bit deterministic\n");
        goto cleanup;
    }

    fdtd_sim_config_t invalid = cfg;
    invalid.n_steps = 0;

    if (fdtd_simulate_echo(&invalid, a, NULL) == 0) {
        fprintf(stderr, "FAIL: simulation accepted zero steps\n");
        goto cleanup;
    }

    invalid = cfg;
    invalid.probe_i = cfg.nx + 1;

    if (fdtd_simulate_echo(&invalid, a, NULL) == 0) {
        fprintf(stderr, "FAIL: simulation accepted an out-of-bounds probe\n");
        goto cleanup;
    }

    if (fdtd_simulate_echo(NULL, a, NULL) == 0) {
        fprintf(stderr, "FAIL: simulation accepted a NULL config\n");
        goto cleanup;
    }

    if (fdtd_simulate_echo(&cfg, NULL, NULL) == 0) {
        fprintf(stderr, "FAIL: simulation accepted a NULL output buffer\n");
        goto cleanup;
    }

    printf("simulate_unit: PASS\n");
    printf("samples: %zu\n", cfg.n_steps);
    printf("peak: %.9e\n", peak);
    printf("early peak: %.9e\n", early_peak);
    printf("sample rate: %.9e Hz\n", fs_a);
    printf("determinism: bit-for-bit identical\n");

    result = 0;

cleanup:
    free(a);
    free(b);
    return result;
}