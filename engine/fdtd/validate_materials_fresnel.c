#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "cpml.h"
#include "materials.h"
#include "source.h"
#include "probe.h"

int main(void) {
    const size_t NX = 210, NY = 160, NZ = 160;
    const double dx = 1e-3;
    const size_t NPML = 8;
    const size_t n_steps = 420;

    size_t i_src = 20;
    size_t j_center = NY / 2, k_center = NZ / 2;
    size_t i_probe = 100;       /* 80 cells from source, ~4 wavelengths */
    size_t i_interface = 150;   /* 130 cells from source, ~6.5 wavelengths */
    size_t i_mat_end = NX - NPML;

    /* small aperture, sized so its own Fraunhofer distance (2*D^2/lambda)
     * stays well under the 80-cell probe distance -- unlike the 65-wide
     * aperture that broke the last test */
    const int APERTURE_HALF = 13; /* width 27, R_ff ~= 73 cells < 80 */
    size_t n_side = (size_t)(2 * APERTURE_HALF + 1);
    size_t n_sources = n_side * n_side;

    double eps_r = 4.0;
    double sigma = 0.0;

    double n1 = 1.0, n2 = sqrt(eps_r);
    double R_theory = (n1 - n2) / (n1 + n2);
    double R_theory_db = 20.0 * log10(fabs(R_theory));

    fdtd_grid_t *g_vac = fdtd_grid_create(NX, NY, NZ, dx, 0.99);
    fdtd_grid_t *g_mat = fdtd_grid_create(NX, NY, NZ, dx, 0.99);
    if (!g_vac || !g_mat) { fprintf(stderr, "grid alloc failed\n"); return 1; }

    cpml_t *pml_vac = cpml_create(g_vac, cpml_default_params(NPML));
    cpml_t *pml_mat = cpml_create(g_mat, cpml_default_params(NPML));
    if (!pml_vac || !pml_mat) { fprintf(stderr, "cpml_create failed\n"); return 1; }

    materials_t *mat = materials_create(g_mat);
    if (!mat) { fprintf(stderr, "materials_create failed\n"); return 1; }

    if (materials_set_box(mat, g_mat, i_interface, i_mat_end + 1, 0, NY + 1, 0, NZ + 1,
                           eps_r, sigma) != 0) {
        fprintf(stderr, "materials_set_box failed\n");
        return 1;
    }

    double t0 = 40.0 * g_vac->dt;
    double tau = 8.0 * g_vac->dt;
    source_gaussian_deriv_params_t src_params = { t0, tau };

    source_t *aperture = malloc(n_sources * sizeof(source_t));
    if (!aperture) { fprintf(stderr, "malloc failed\n"); return 1; }

    double per_source_amplitude = 1.0 / (double)n_sources;

    size_t s = 0;
    for (int dj = -APERTURE_HALF; dj <= APERTURE_HALF; dj++) {
        for (int dk = -APERTURE_HALF; dk <= APERTURE_HALF; dk++) {
            size_t j = j_center + (size_t)dj;
            size_t k = k_center + (size_t)dk;
            aperture[s] = source_make_point(i_src, j, k, FIELD_EZ, SRC_SOFT,
                                             per_source_amplitude,
                                             source_waveform_gaussian_derivative,
                                             &src_params);
            s++;
        }
    }

    probe_td_t probe_vac, probe_mat;
    if (probe_td_create(&probe_vac, i_probe, j_center, k_center, FIELD_EZ, n_steps) != 0 ||
        probe_td_create(&probe_mat, i_probe, j_center, k_center, FIELD_EZ, n_steps) != 0) {
        fprintf(stderr, "probe_td_create failed\n");
        return 1;
    }

    for (size_t n = 0; n < n_steps; n++) {
        double t = (double)n * g_vac->dt;

        fdtd_update_h_cpml(g_vac, pml_vac);
        fdtd_update_h_cpml(g_mat, pml_mat);

        for (size_t si = 0; si < n_sources; si++) {
            source_apply(&aperture[si], g_vac, t);
            source_apply(&aperture[si], g_mat, t);
        }

        fdtd_update_e_cpml(g_vac, pml_vac);
        fdtd_update_e_cpml_materials(g_mat, pml_mat, mat);

        probe_td_record(&probe_vac, g_vac);
        probe_td_record(&probe_mat, g_mat);

        if (n % 50 == 0) {
            fprintf(stderr, "step %zu / %zu\n", n, n_steps);
        }
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

    size_t incident_start = 150, incident_end = 230;
    size_t reflect_start = 300, reflect_end = 419;

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

    printf("materials/Fresnel validation -- large domain (eps_r=%.2f, sigma=%.2f)\n", eps_r, sigma);
    printf("  grid: %zux%zux%zu, aperture: %zux%zu (~1.35 wavelengths wide)\n", NX, NY, NZ, n_side, n_side);
    printf("  source-probe distance: %zu cells, source-interface distance: %zu cells\n",
           i_probe - i_src, i_interface - i_src);
    printf("  incident peak  (steps %zu-%zu): %.6e\n", incident_start, incident_end, incident_peak);
    printf("  reflected peak (steps %zu-%zu): %.6e\n", reflect_start, reflect_end, reflected_peak);
    printf("  measured |R| = %.6f  (%.3f dB)\n", R_measured, R_measured_db);
    printf("  theory   |R| = %.6f  (%.3f dB)\n", fabs(R_theory), R_theory_db);
    printf("  error: %.3f dB, %.2f%% relative to theoretical |R|\n",
           R_measured_db - R_theory_db,
           100.0 * fabs(R_measured - fabs(R_theory)) / fabs(R_theory));

    free(diff);
    free(aperture);
    probe_td_destroy(&probe_vac);
    probe_td_destroy(&probe_mat);
    materials_destroy(mat);
    cpml_destroy(pml_vac);
    cpml_destroy(pml_mat);
    fdtd_grid_destroy(g_vac);
    fdtd_grid_destroy(g_mat);

    return 0;
}