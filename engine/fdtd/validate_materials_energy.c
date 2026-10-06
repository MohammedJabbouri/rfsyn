#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "cpml.h"
#include "materials.h"
#include "source.h"

static void snapshot_hy_hz(const fdtd_grid_t *g, size_t i_plane,
                            double *hy_buf, double *hz_buf)
{
    size_t nzp1 = g->nz + 1;
    for (size_t j = 0; j <= g->ny; j++) {
        for (size_t k = 0; k <= g->nz; k++) {
            size_t idx = fdtd_index(g, i_plane, j, k);
            hy_buf[j * nzp1 + k] = g->hy[idx];
            hz_buf[j * nzp1 + k] = g->hz[idx];
        }
    }
}

static double flux_power(const fdtd_grid_t *g, size_t i_plane,
                          const double *hy_prev, const double *hz_prev,
                          double dx)
{
    double power = 0.0;
    size_t nzp1 = g->nz + 1;
    for (size_t j = 0; j <= g->ny; j++) {
        for (size_t k = 0; k <= g->nz; k++) {
            size_t idx = fdtd_index(g, i_plane, j, k);
            double hy_avg = 0.5 * (hy_prev[j * nzp1 + k] + g->hy[idx]);
            double hz_avg = 0.5 * (hz_prev[j * nzp1 + k] + g->hz[idx]);
            double sx = g->ey[idx] * hz_avg - g->ez[idx] * hy_avg;
            power += sx;
        }
    }
    return power * dx * dx;
}

static double flux_power_reflected(const fdtd_grid_t *g_mat, const fdtd_grid_t *g_vac,
                                    size_t i_plane,
                                    const double *hy_prev_mat, const double *hz_prev_mat,
                                    const double *hy_prev_vac, const double *hz_prev_vac,
                                    double dx)
{
    double power = 0.0;
    size_t nzp1 = g_mat->nz + 1;
    for (size_t j = 0; j <= g_mat->ny; j++) {
        for (size_t k = 0; k <= g_mat->nz; k++) {
            size_t idx_mat = fdtd_index(g_mat, i_plane, j, k);
            size_t idx_vac = fdtd_index(g_vac, i_plane, j, k);

            double hy_avg_mat = 0.5 * (hy_prev_mat[j * nzp1 + k] + g_mat->hy[idx_mat]);
            double hz_avg_mat = 0.5 * (hz_prev_mat[j * nzp1 + k] + g_mat->hz[idx_mat]);
            double hy_avg_vac = 0.5 * (hy_prev_vac[j * nzp1 + k] + g_vac->hy[idx_vac]);
            double hz_avg_vac = 0.5 * (hz_prev_vac[j * nzp1 + k] + g_vac->hz[idx_vac]);

            double hy_refl = hy_avg_mat - hy_avg_vac;
            double hz_refl = hz_avg_mat - hz_avg_vac;
            double ey_refl = g_mat->ey[idx_mat] - g_vac->ey[idx_vac];
            double ez_refl = g_mat->ez[idx_mat] - g_vac->ez[idx_vac];

            double sx = ey_refl * hz_refl - ez_refl * hy_refl;
            power += sx;
        }
    }
    return power * dx * dx;
}

int main(void) {
    const size_t N = 110;
    const double dx = 1e-3;
    const size_t NPML = 8;
    const size_t n_steps = 400;

    size_t i_src = 30, j_src = 55, k_src = 55;

    size_t i_plane1 = 48;
    size_t i_interface = 50;
    size_t i_plane2 = 52;
    size_t i_mat_end = N - NPML;

    double eps_r = 4.0;
    double sigma = 0.0;

    fdtd_grid_t *g_vac = fdtd_grid_create(N, N, N, dx, 0.99);
    fdtd_grid_t *g_mat = fdtd_grid_create(N, N, N, dx, 0.99);
    if (!g_vac || !g_mat) { fprintf(stderr, "grid alloc failed\n"); return 1; }

    cpml_t *pml_vac = cpml_create(g_vac, cpml_default_params(NPML));
    cpml_t *pml_mat = cpml_create(g_mat, cpml_default_params(NPML));
    if (!pml_vac || !pml_mat) { fprintf(stderr, "cpml_create failed\n"); return 1; }

    materials_t *mat = materials_create(g_mat);
    if (!mat) { fprintf(stderr, "materials_create failed\n"); return 1; }

    if (materials_set_box(mat, g_mat, i_interface, i_mat_end + 1, 0, N + 1, 0, N + 1, eps_r, sigma) != 0) {
        fprintf(stderr, "materials_set_box failed\n");
        return 1;
    }

    double t0 = 40.0 * g_vac->dt;
    double tau = 8.0 * g_vac->dt;
    source_gaussian_deriv_params_t src_params = { t0, tau };
    source_t src = source_make_point(i_src, j_src, k_src, FIELD_EZ, SRC_SOFT, 1.0,
                                      source_waveform_gaussian_derivative, &src_params);

    size_t plane_size = (N + 1) * (N + 1);
    double *hy_prev_vac_p1 = malloc(plane_size * sizeof(double));
    double *hz_prev_vac_p1 = malloc(plane_size * sizeof(double));
    double *hy_prev_mat_p1 = malloc(plane_size * sizeof(double));
    double *hz_prev_mat_p1 = malloc(plane_size * sizeof(double));
    double *hy_prev_mat_p2 = malloc(plane_size * sizeof(double));
    double *hz_prev_mat_p2 = malloc(plane_size * sizeof(double));

    if (!hy_prev_vac_p1 || !hz_prev_vac_p1 || !hy_prev_mat_p1 ||
        !hz_prev_mat_p1 || !hy_prev_mat_p2 || !hz_prev_mat_p2) {
        fprintf(stderr, "malloc failed\n");
        return 1;
    }

    double energy_incident = 0.0, energy_reflected = 0.0, energy_transmitted = 0.0;

    FILE *out = fopen("materials_energy_validation.csv", "w");
    fprintf(out, "step,t,energy_incident_cum,energy_reflected_cum,energy_transmitted_cum\n");

    for (size_t n = 0; n < n_steps; n++) {
        double t = (double)n * g_vac->dt;

        snapshot_hy_hz(g_vac, i_plane1, hy_prev_vac_p1, hz_prev_vac_p1);
        snapshot_hy_hz(g_mat, i_plane1, hy_prev_mat_p1, hz_prev_mat_p1);
        snapshot_hy_hz(g_mat, i_plane2, hy_prev_mat_p2, hz_prev_mat_p2);

        fdtd_update_h_cpml(g_vac, pml_vac);
        fdtd_update_h_cpml(g_mat, pml_mat);

        double p_inc   = flux_power(g_vac, i_plane1, hy_prev_vac_p1, hz_prev_vac_p1, dx);
        double p_refl  = flux_power_reflected(g_mat, g_vac, i_plane1,
                                               hy_prev_mat_p1, hz_prev_mat_p1,
                                               hy_prev_vac_p1, hz_prev_vac_p1, dx);
        double p_trans = flux_power(g_mat, i_plane2, hy_prev_mat_p2, hz_prev_mat_p2, dx);

        energy_incident    += p_inc * g_vac->dt;
        energy_reflected   += p_refl * g_vac->dt;
        energy_transmitted += p_trans * g_vac->dt;

        source_apply(&src, g_vac, t);
        source_apply(&src, g_mat, t);

        fdtd_update_e_cpml(g_vac, pml_vac);
        fdtd_update_e_cpml_materials(g_mat, pml_mat, mat);

        fprintf(out, "%zu,%.9e,%.9e,%.9e,%.9e\n", n, t, energy_incident, energy_reflected, energy_transmitted);

        if (n % 50 == 0) fprintf(stderr, "step %zu / %zu\n", n, n_steps);
    }

    fclose(out);

    double refl_mag = fabs(energy_reflected);
    double sum = refl_mag + energy_transmitted;
    double mismatch_pct = 100.0 * fabs(sum - energy_incident) / energy_incident;

    printf("materials energy-conservation validation (eps_r=%.2f, sigma=%.2f)\n", eps_r, sigma);
    printf("  incident energy    : %.6e\n", energy_incident);
    printf("  reflected energy    : %.6e (raw signed: %.6e)\n", refl_mag, energy_reflected);
    printf("  transmitted energy  : %.6e\n", energy_transmitted);
    printf("  reflected+transmitted: %.6e\n", sum);
    printf("  mismatch vs incident: %.3f%%\n", mismatch_pct);

    free(hy_prev_vac_p1); free(hz_prev_vac_p1);
    free(hy_prev_mat_p1); free(hz_prev_mat_p1);
    free(hy_prev_mat_p2); free(hz_prev_mat_p2);
    materials_destroy(mat);
    cpml_destroy(pml_vac);
    cpml_destroy(pml_mat);
    fdtd_grid_destroy(g_vac);
    fdtd_grid_destroy(g_mat);

    return 0;
}