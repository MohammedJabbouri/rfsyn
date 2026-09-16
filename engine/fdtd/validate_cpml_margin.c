#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "cpml.h"

static double source_waveform(double t, double t0, double tau) {
    double arg = (t - t0) / tau;
    return -2.0 * arg * exp(-arg * arg);
}

#define N_PROBES 5

int main(void) {
    const size_t N_TEST = 50;
    const size_t N_REF  = 150;

    const double dx = 1e-3;
    const size_t NPML = 8;
    const size_t n_steps = 200;

    const size_t offsets[N_PROBES] = {6, 9, 12, 15, 16};

    size_t center_test = N_TEST / 2;
    size_t center_ref  = N_REF / 2;
    size_t pml_interface_test = N_TEST - NPML;

    fdtd_grid_t *g_ref  = fdtd_grid_create(N_REF, N_REF, N_REF, dx, 0.99);
    fdtd_grid_t *g_test = fdtd_grid_create(N_TEST, N_TEST, N_TEST, dx, 0.99);
    if (!g_ref || !g_test) { fprintf(stderr, "grid alloc failed\n"); return 1; }

    if (fabs(g_ref->dt - g_test->dt) > 1e-30) {
        fprintf(stderr, "dt mismatch -- unexpected\n");
        return 1;
    }

    cpml_t *pml = cpml_create(g_test, cpml_default_params(NPML));
    if (!pml) { fprintf(stderr, "cpml_create failed\n"); return 1; }

    double t0 = 40.0 * g_ref->dt;
    double tau = 8.0 * g_ref->dt;

    size_t probe_test_i[N_PROBES], probe_ref_i[N_PROBES], dist_to_pml[N_PROBES];
    for (int p = 0; p < N_PROBES; p++) {
        probe_test_i[p] = center_test + offsets[p];
        probe_ref_i[p]  = center_ref + offsets[p];
        dist_to_pml[p]  = pml_interface_test - probe_test_i[p];
    }

    double max_ref_mag[N_PROBES] = {0};
    double sum_db[N_PROBES], max_db[N_PROBES];
    long   count_db[N_PROBES];
    for (int p = 0; p < N_PROBES; p++) { sum_db[p]=0.0; max_db[p]=-1e300; count_db[p]=0; }

    FILE *out = fopen("cpml_margin_validation.csv", "w");
    fprintf(out, "step,t,probe_id,dist_to_pml,probe_ref,probe_test,abs_err,db_err\n");

    const size_t SUMMARY_START = 60;

    for (size_t n = 0; n < n_steps; n++) {
        double t = (double)n * g_ref->dt;
        double src_val = source_waveform(t, t0, tau);

        fdtd_update_h(g_ref);
        fdtd_update_h_cpml(g_test, pml);

        g_ref->ez[fdtd_index(g_ref, center_ref, center_ref, center_ref)] += src_val;
        g_test->ez[fdtd_index(g_test, center_test, center_test, center_test)] += src_val;

        fdtd_update_e(g_ref);
        fdtd_update_e_cpml(g_test, pml);

        for (int p = 0; p < N_PROBES; p++) {
            double p_ref  = g_ref->ez[fdtd_index(g_ref, probe_ref_i[p], center_ref, center_ref)];
            double p_test = g_test->ez[fdtd_index(g_test, probe_test_i[p], center_test, center_test)];

            if (fabs(p_ref) > max_ref_mag[p]) max_ref_mag[p] = fabs(p_ref);

            double abs_err = fabs(p_test - p_ref);
            double db_err = (max_ref_mag[p] > 0.0)
            ? 20.0 * log10(abs_err / max_ref_mag[p] + 1e-300)
            : 0.0;

            fprintf(out, "%zu,%.9e,%d,%zu,%.9e,%.9e,%.9e,%.6f\n", n, t, p, dist_to_pml[p], p_ref, p_test, abs_err, db_err);

            if (n >= SUMMARY_START) {
                sum_db[p] += db_err;
                if (db_err > max_db[p]) max_db[p] = db_err;
                count_db[p]++;
            }
        }
    }
    fclose(out);

    printf("wrote cpml_margin_validation.csv\n\n");
    printf("reflection floor vs. distance from PML interface (steps %zu-%zu):\n", SUMMARY_START, n_steps);
    printf("%-14s %-12s %-12s\n", "dist_to_pml", "mean_db", "worst_db");
    for (int p = 0; p < N_PROBES; p++) {
        double mean_db = sum_db[p] / (double)count_db[p];
        printf("%-14zu %-12.2f %-12.2f\n", dist_to_pml[p], mean_db, max_db[p]);
    }

    cpml_destroy(pml);
    fdtd_grid_destroy(g_test);
    fdtd_grid_destroy(g_ref);
    return 0;
}