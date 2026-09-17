#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "cpml.h"
#include "source.h"
#include "probe.h"

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
        fprintf(stderr, "dt mismatch | unexpected\n");
        return 1;
    }

    cpml_t *pml = cpml_create(g_test, cpml_default_params(NPML));
    if (!pml) { fprintf(stderr, "cpml_create failed\n"); return 1; }

    double t0 = 40.0 * g_ref->dt;
    double tau = 8.0 * g_ref->dt;
    source_gaussian_deriv_params_t src_params = { t0, tau };

    source_t src_ref  = source_make_point(center_ref, center_ref, center_ref,
                                           FIELD_EZ, SRC_SOFT, 1.0,
                                           source_waveform_gaussian_derivative, &src_params);
    source_t src_test = source_make_point(center_test, center_test, center_test,
                                           FIELD_EZ, SRC_SOFT, 1.0,
                                           source_waveform_gaussian_derivative, &src_params);

    double freq_hz = 1.5e10;

    probe_dft_t probes_ref[N_PROBES];
    probe_dft_t probes_test[N_PROBES];
    size_t dist_to_pml[N_PROBES];

    for (size_t p = 0; p < N_PROBES; p++) {
        size_t i_ref  = center_ref + offsets[p];
        size_t i_test = center_test + offsets[p];
        dist_to_pml[p] = pml_interface_test - i_test;

        if (probe_dft_create(&probes_ref[p], i_ref, center_ref, center_ref,
                              FIELD_EZ, &freq_hz, 1) != 0) {
            fprintf(stderr, "probe_dft_create (ref) failed\n");
            return 1;
        }
        if (probe_dft_create(&probes_test[p], i_test, center_test, center_test,
                              FIELD_EZ, &freq_hz, 1) != 0) {
            fprintf(stderr, "probe_dft_create (test) failed\n");
            return 1;
        }
    }

    for (size_t n = 0; n < n_steps; n++) {
        double t = (double)n * g_ref->dt;

        fdtd_update_h(g_ref);
        fdtd_update_h_cpml(g_test, pml);

        source_apply(&src_ref, g_ref, t);
        source_apply(&src_test, g_test, t);

        fdtd_update_e(g_ref);
        fdtd_update_e_cpml(g_test, pml);

        for (size_t p = 0; p < N_PROBES; p++) {
            probe_dft_record(&probes_ref[p], g_ref, t, g_ref->dt);
            probe_dft_record(&probes_test[p], g_test, t, g_test->dt);
        }
    }

    printf("cross-check: source.c/probe.c vs. previously-validated margin result\n");
    printf("frequency = %.3e Hz\n\n", freq_hz);
    printf("%-14s %-14s %-14s %-12s\n", "dist_to_pml", "mag_ref", "diff_mag", "error_db");

    for (size_t p = 0; p < N_PROBES; p++) {
        double mag_ref, phase_ref;
        probe_dft_result(&probes_ref[p], 0, &mag_ref, &phase_ref);

        double diff_re = probes_test[p].acc_re[0] - probes_ref[p].acc_re[0];
        double diff_im = probes_test[p].acc_im[0] - probes_ref[p].acc_im[0];
        double diff_mag = sqrt(diff_re * diff_re + diff_im * diff_im);

        double error_db = (mag_ref > 0.0)
            ? 20.0 * log10(diff_mag / mag_ref + 1e-300)
            : 0.0;

        printf("%-14zu %-14.6e %-14.6e %-12.3f\n",
               dist_to_pml[p], mag_ref, diff_mag, error_db);
    }

    for (size_t p = 0; p < N_PROBES; p++) {
        probe_dft_destroy(&probes_ref[p]);
        probe_dft_destroy(&probes_test[p]);
    }
    cpml_destroy(pml);
    fdtd_grid_destroy(g_test);
    fdtd_grid_destroy(g_ref);

    return 0;
}