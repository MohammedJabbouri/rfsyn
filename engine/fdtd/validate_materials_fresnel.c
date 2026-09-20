#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "cpml.h"
#include "materials.h"
#include "source.h"
#include "probe.h"

int main(void) {
    const size_t N = 110;
    const double dx = 1e-3;
    const size_t NPML = 8;
    const size_t n_steps = 260;

    size_t i_src = 20, j_src = 55, k_src = 55;
    size_t i_probe = 70;
    size_t i_interface = 90;
    size_t i_mat_end = N - NPML;

    double eps_r = 4.0;
    double sigma = 0.0;

    double n1 = 1.0, n2 = sqrt(eps_r);
    double R_theory = (n1 - n2) / (n1 + n2);
    double R_theory_db = 20.0 * log10(fabs(R_theory));

    fdtd_grid_t *g_vac = fdtd_grid_create(N, N, N, dx, 0.99);
    fdtd_grid_t *g_mat = fdtd_grid_create(N, N, N, dx, 0.99);
    if (!g_vac || !g_mat) { fprintf(stderr, "grid alloc failed\n"); return 1; }

    cpml_t *pml_vac = cpml_create(g_vac, cpml_default_params(NPML));
    cpml_t *pml_mat = cpml_create(g_mat, cpml_default_params(NPML));
    if (!pml_vac || !pml_mat) { fprintf(stderr, "cpml_create failed\n"); return 1; }

    materials_t *mat = materials_create(g_mat);
    if (!mat) { fprintf(stderr, "materials_create failed\n"); return 1; }

    if (materials_set_box(mat, g_mat, i_interface, i_mat_end + 1, 0, N + 1, 0, N + 1,
                           eps_r, sigma) != 0) {
        fprintf(stderr, "materials_set_box failed\n");
        return 1;
    }

    double t0 = 40.0 * g_vac->dt;
    double tau = 8.0 * g_vac->dt;
    source_gaussian_deriv_params_t src_params = { t0, tau };
    source_t src = source_make_point(i_src, j_src, k_src, FIELD_EZ, SRC_SOFT, 1.0, source_waveform_gaussian_derivative, &src_params);

    probe_td_t probe_vac, probe_mat;
    if (probe_td_create(&probe_vac, i_probe, j_src, k_src, FIELD_EZ, n_steps) != 0 ||
        probe_td_create(&probe_mat, i_probe, j_src, k_src, FIELD_EZ, n_steps) != 0) {
        fprintf(stderr, "probe_td_create failed\n");
        return 1;
    }

    for (size_t n = 0; n < n_steps; n++) {
        double t = (double)n * g_vac->dt;

        fdtd_update_h_cpml(g_vac, pml_vac);
        fdtd_update_h_cpml(g_mat, pml_mat);

        source_apply(&src, g_vac, t);
        source_apply(&src, g_mat, t);

        fdtd_update_e_cpml(g_vac, pml_vac);
        fdtd_update_e_cpml_materials(g_mat, pml_mat, mat);

        probe_td_record(&probe_vac, g_vac);
        probe_td_record(&probe_mat, g_mat);
    }

    double *diff = malloc(n_steps * sizeof(double));
    if (!diff) { fprintf(stderr, "malloc failed\n"); return 1; }

    FILE *out = fopen("materials_fresnel_validation.csv", "w");
    fprintf(out, "step,t,probe_vacuum,probe_material,diff\n");
    for (size_t n = 0; n < n_steps; n++) {
        double t = (double)n * g_vac->dt;
        diff[n] = probe_mat.buffer[n] - probe_vac.buffer[n];
        fprintf(out, "%zu,%.9e,%.9e,%.9e,%.9e\n",
                n, t, probe_vac.buffer[n], probe_mat.buffer[n], diff[n]);
    }
    fclose(out);

    size_t incident_start = 90, incident_end = 180;

    size_t reflect_start = 181, reflect_end = 259;

    if (incident_end >= n_steps) incident_end = n_steps - 1;
    if (reflect_end >= n_steps) reflect_end = n_steps - 1;

    double incident_peak = 0.0;
    for (size_t n = incident_start; n <= incident_end; n++) {
        double v = fabs(probe_vac.buffer[n]);
        if (v > incident_peak) incident_peak = v;
    }

    double reflected_peak = 0.0;
    for (size_t n = reflect_start; n <= reflect_end; n++) {
        double v = fabs(diff[n]);
        if (v > reflected_peak) reflected_peak = v;
    }

    double R_measured = reflected_peak / incident_peak;
    double R_measured_db = 20.0 * log10(R_measured);

    printf("materials/Fresnel validation (eps_r=%.2f, sigma=%.2f)\n", eps_r, sigma);
    printf("  source-probe distance: %zu cells, source-interface distance: %zu cells\n",
           i_probe - i_src, i_interface - i_src);
    printf("  incident peak  (steps %zu-%zu): %.6e\n", incident_start, incident_end, incident_peak);
    printf("  reflected peak (steps %zu-%zu): %.6e\n", reflect_start, reflect_end, reflected_peak);
    printf("  measured |R| = %.6f  (%.3f dB)\n", R_measured, R_measured_db);
    printf("  theoritical   |R| = %.6f  (%.3f dB)\n", fabs(R_theory), R_theory_db);
    printf("  error: %.3f dB, %.2f%% relative to theoretical |R|\n", R_measured_db - R_theory_db, 100.0 * fabs(R_measured - fabs(R_theory)) / fabs(R_theory));

    free(diff);
    probe_td_destroy(&probe_vac);
    probe_td_destroy(&probe_mat);
    materials_destroy(mat);
    cpml_destroy(pml_vac);
    cpml_destroy(pml_mat);
    fdtd_grid_destroy(g_vac);
    fdtd_grid_destroy(g_mat);

    return 0;
}