#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "cpml.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double source_waveform(double t, double t0, double tau) {
    double arg = (t - t0) / tau;
    return -2.0 * arg * exp(-arg * arg);
}

static void inject_plane_wave(fdtd_grid_t *g, size_t i_src, double val) {
    for (size_t j = 1; j < g->ny; j++) {
        for (size_t k = 0; k < g->nz; k++) {
            g->ez[fdtd_index(g, i_src, j, k)] += val;
        }
    }
}

int main(void) {
    const double dx = 1e-3;
    const size_t NPML = 8;

    const size_t NX_TEST = 60, NY_TEST = 20, NZ_TEST = 20;
    size_t src_test = 15, probe_test_i = 45;

    const size_t NX_REF = 400, NY_REF = 20, NZ_REF = 20;
    size_t src_ref = 200, probe_ref_i = 230;

    const size_t n_steps = 250;

    fdtd_grid_t *g_ref = fdtd_grid_create(NX_REF, NY_REF, NZ_REF, dx, 0.99);
    fdtd_grid_t *g_test = fdtd_grid_create(NX_TEST, NY_TEST, NZ_TEST, dx, 0.99);
    if (!g_ref || !g_test) { fprintf(stderr, "grid alloc failed\n"); return 1; }

    if (fabs(g_ref->dt - g_test->dt) > 1e-30) {
        fprintf(stderr, "dt mismatch -- unexpected\n");
        return 1;
    }

    cpml_t *pml = cpml_create(g_test, cpml_default_params(NPML));
    if (!pml) { fprintf(stderr, "cpml_create failed\n"); return 1; }

    double t0 = 40.0 * g_ref->dt;
    double tau = 8.0 * g_ref->dt;

    size_t j_probe_test = NY_TEST / 2, k_probe_test = NZ_TEST / 2;
    size_t j_probe_ref  = NY_REF / 2,  k_probe_ref  = NZ_REF / 2;

    FILE *out = fopen("cpml_planewave_validation.csv", "w");
    fprintf(out, "step, t, probe_ref,probe_test,abs_err,db_err\n");

    double max_ref_mag = 0.0;

    for (size_t n = 0; n < n_steps; n++) {
        double t = (double)n * g_ref->dt;
        double src_val = source_waveform(t, t0, tau);

        fdtd_update_h(g_ref);
        fdtd_update_h_cpml(g_test, pml);

        inject_plane_wave(g_ref, src_ref, src_val);
        inject_plane_wave(g_test, src_test, src_val);

        fdtd_update_e(g_ref);
        fdtd_update_e_cpml(g_test, pml);

        double p_ref  = g_ref->ez[fdtd_index(g_ref, probe_ref_i, j_probe_ref, k_probe_ref)];
        double p_test = g_test->ez[fdtd_index(g_test, probe_test_i, j_probe_test, k_probe_test)];

        if (fabs(p_ref) > max_ref_mag) max_ref_mag = fabs(p_ref);

        double abs_err = fabs(p_test - p_ref);
        double db_err = (max_ref_mag > 0.0)
        ? 20.0 * log10(abs_err / max_ref_mag + 1e-300)
        : 0.0;

        fprintf(out, "%zu,%.9e,%.9e,%.9e,%.9e,%.6f\n", n, t, p_ref, p_test, abs_err, db_err);
    }

    fclose(out);
    cpml_destroy(pml);
    fdtd_grid_destroy(g_test);
    fdtd_grid_destroy(g_ref);

    printf("wrote cpml_planewave_validation.csv\n");
    printf("normal-incidence test -- expect floor near r0 target (1e-6 -> -120 dB)\n");
    return 0;
}