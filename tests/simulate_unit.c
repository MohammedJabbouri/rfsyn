#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../engine/fdtd/simulate.h"

int main(void) {
    fdtd_sim_config_t cfg = {
        .nx = 50, .ny = 50, .nz = 50,
        .dx = 1e-3, .courant_safety = 0.99,
        .n_steps = 200, .npml = 8,
        .src_i = 25, .src_j = 25, .src_k = 25,
        .src_component = FIELD_EZ,
        .t0_over_dt = 40.0, .tau_over_dt = 8.0,
        .probe_i = 31, .probe_j = 25, .probe_k = 25,
        .probe_component = FIELD_EZ,
    };

    double *a = malloc(cfg.n_steps * sizeof(double));
    double *b = malloc(cfg.n_steps * sizeof(double));
    if (!a || !b) { fprintf(stderr, "malloc failed\n"); return 1; }

    double fs = 0.0;
    assert(fdtd_simulate_echo(&cfg, a, &fs) == 0);
    assert(fdtd_simulate_echo(&cfg, b, NULL) == 0);

    assert(memcmp(a, b, cfg.n_steps * sizeof(double)) == 0);

    double peak = 0.0, early = 0.0;
    for (size_t i = 0; i < cfg.n_steps; i++) {
        assert(isfinite(a[i]));
        double v = fabs(a[i]);
        if (v > peak) peak = v;
        if (i < 20 && v > early) early = v;
    }
    assert(peak > 1e-9);
    assert(early < peak * 1e-3);
    assert(isfinite(fs) && fs > 0.0);

    printf("simulate_unit: PASS (peak=%.3e, fs=%.3e Hz)\n", peak, fs);
    free(a); free(b);
    return 0;
}