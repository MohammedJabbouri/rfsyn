#include "cpml.h"
#include "../core/constants.h"
#include <stdlib.h>
#include <math.h>

cpml_params_t cpml_default_params(size_t npml) {
    cpml_params_t p;
    p.npml = npml;
    p.m = 3.0;
    p.ma = 1.0;
    p.kappa_max = 7.0;
    p.alpha_max = 0.05;
    p.r0 = 1e-6;
    return p;
}

static int build_axis_profiles(size_t n, size_t npml, double dx, double m, double ma,
    double kappa_max, double alpha_max, double sigma_max, double dt,
    double **kinv_m_out, double **b_m_out, double **a_m_out,
    double **kinv_h_out, double **b_h_out, double **a_h_out)
{
    size_t nm = n + 1;
    size_t nh = n;

    double *KM = malloc(nm * sizeof(double));
    double *BM = malloc(nm * sizeof(double));
    double *AM = malloc(nm * sizeof(double));
    double *KH = malloc(nh * sizeof(double));
    double *BH = malloc(nh * sizeof(double));
    double *AH = malloc(nh * sizeof(double));

    if (!KM || !BM || !AM || !KH || !BH || !AH) {
        free(KM); free(BM); free(AM); free(KH); free(BH); free(AH);
        return -1;
    }

    double d = (double)npml * dx;
    double domain_len = (double)n * dx;

    for (size_t idx = 0; idx < nm; idx++) {
        double x = (double)idx * dx;
        double rho = -1.0;
        if (x < d) rho = d - x;
        else if (x > domain_len - d) rho = x - (domain_len - d);

        if (rho < 0.0) {
            KM[idx] = 1.0; BM[idx] = 1.0; AM[idx] = 0.0;
        } else {
            double r = rho / d; if (r > 1.0) r = 1.0;
            double sigma = sigma_max * pow(r, m);
            double kappa = 1.0 + (kappa_max - 1.0) * pow(r, m);
            double alpha = alpha_max * pow(1.0 - r, ma);
            double b = exp(-(sigma / kappa + alpha) * dt / EPS0);
            double a = (sigma > 0.0)
                ? sigma * (b - 1.0) / (dx * kappa * (sigma + kappa * alpha))
                : 0.0;
            KM[idx] = 1.0 / kappa; BM[idx] = b; AM[idx] = a;
        }
    }

    for (size_t idx = 0; idx < nh; idx++) {
        double x = ((double)idx + 0.5) * dx;
        double rho = -1.0;
        if (x < d) rho = d - x;
        else if (x > domain_len - d) rho = x - (domain_len - d);

        if (rho < 0.0) {
            KH[idx] = 1.0; BH[idx] = 1.0; AH[idx] = 0.0;
        } else {
            double r = rho / d; if (r > 1.0) r = 1.0;
            double sigma = sigma_max * pow(r, m);
            double kappa = 1.0 + (kappa_max - 1.0) * pow(r, m);
            double alpha = alpha_max * pow(1.0 - r, ma);
            double b = exp(-(sigma / kappa + alpha) * dt / EPS0);
            double a = (sigma > 0.0)
                ? sigma * (b - 1.0) / (dx * kappa * (sigma + kappa * alpha))
                : 0.0;
            KH[idx] = 1.0 / kappa; BH[idx] = b; AH[idx] = a;
        }
    }

    *kinv_m_out = KM; *b_m_out = BM; *a_m_out = AM;
    *kinv_h_out = KH; *b_h_out = BH; *a_h_out = AH;
    return 0;
}

cpml_t *cpml_create(const fdtd_grid_t *g, cpml_params_t params) {
    if (!g || params.npml == 0) return NULL;
    if (2 * params.npml >= g->nx || 2 * params.npml >= g->ny || 2 * params.npml >= g->nz) {
        return NULL;
    }

    cpml_t *pml = calloc(1, sizeof(cpml_t));
    if (!pml) return NULL;

    pml->npml = params.npml;
    pml->dx = g->dx;

    double eta0 = sqrt(MU0 / EPS0);
    double sigma_max = -(params.m + 1.0) * log(params.r0) /
                        (2.0 * eta0 * (double)params.npml * g->dx);

    int ok = 1;
    ok &= (build_axis_profiles(g->nx, params.npml, g->dx, params.m, params.ma, params.kappa_max, params.alpha_max, sigma_max, g->dt, &pml->kappa_inv_m_x, &pml->b_m_x, &pml->a_m_x, &pml->kappa_inv_h_x, &pml->b_h_x, &pml->a_h_x) == 0);
    ok &= (build_axis_profiles(g->ny, params.npml, g->dx, params.m, params.ma, params.kappa_max, params.alpha_max, sigma_max, g->dt, &pml->kappa_inv_m_y, &pml->b_m_y, &pml->a_m_y, &pml->kappa_inv_h_y, &pml->b_h_y, &pml->a_h_y) == 0);
    ok &= (build_axis_profiles(g->nz, params.npml, g->dx, params.m, params.ma, params.kappa_max, params.alpha_max, sigma_max, g->dt, &pml->kappa_inv_m_z, &pml->b_m_z, &pml->a_m_z, &pml->kappa_inv_h_z, &pml->b_h_z, &pml->a_h_z) == 0);

    if (!ok) { cpml_destroy(pml); return NULL; }

    size_t total = (g->nx + 1) * (g->ny + 1) * (g->nz + 1);
    pml->psi_hx_y = calloc(total, sizeof(double));
    pml->psi_hx_z = calloc(total, sizeof(double));
    pml->psi_hy_z = calloc(total, sizeof(double));
    pml->psi_hy_x = calloc(total, sizeof(double));
    pml->psi_hz_x = calloc(total, sizeof(double));
    pml->psi_hz_y = calloc(total, sizeof(double));
    pml->psi_ex_y = calloc(total, sizeof(double));
    pml->psi_ex_z = calloc(total, sizeof(double));
    pml->psi_ey_z = calloc(total, sizeof(double));
    pml->psi_ey_x = calloc(total, sizeof(double));
    pml->psi_ez_x = calloc(total, sizeof(double));
    pml->psi_ez_y = calloc(total, sizeof(double));

    if (!pml->psi_hx_y || !pml->psi_hx_z || !pml->psi_hy_z || !pml->psi_hy_x ||
        !pml->psi_hz_x || !pml->psi_hz_y || !pml->psi_ex_y || !pml->psi_ex_z ||
        !pml->psi_ey_z || !pml->psi_ey_x || !pml->psi_ez_x || !pml->psi_ez_y) {
        cpml_destroy(pml);
        return NULL;
    }

    return pml;
}

void cpml_destroy(cpml_t *pml) {
    if (!pml) return;
    free(pml->kappa_inv_m_x); free(pml->b_m_x); free(pml->a_m_x);
    free(pml->kappa_inv_m_y); free(pml->b_m_y); free(pml->a_m_y);
    free(pml->kappa_inv_m_z); free(pml->b_m_z); free(pml->a_m_z);
    free(pml->kappa_inv_h_x); free(pml->b_h_x); free(pml->a_h_x);
    free(pml->kappa_inv_h_y); free(pml->b_h_y); free(pml->a_h_y);
    free(pml->kappa_inv_h_z); free(pml->b_h_z); free(pml->a_h_z);
    free(pml->psi_hx_y); free(pml->psi_hx_z);
    free(pml->psi_hy_z); free(pml->psi_hy_x);
    free(pml->psi_hz_x); free(pml->psi_hz_y);
    free(pml->psi_ex_y); free(pml->psi_ex_z);
    free(pml->psi_ey_z); free(pml->psi_ey_x);
    free(pml->psi_ez_x); free(pml->psi_ez_y);
    free(pml);
}

#define EX(i,j,k) g->ex[fdtd_index(g,(i),(j),(k))]
#define EY(i,j,k) g->ey[fdtd_index(g,(i),(j),(k))]
#define EZ(i,j,k) g->ez[fdtd_index(g,(i),(j),(k))]
#define HX(i,j,k) g->hx[fdtd_index(g,(i),(j),(k))]
#define HY(i,j,k) g->hy[fdtd_index(g,(i),(j),(k))]
#define HZ(i,j,k) g->hz[fdtd_index(g,(i),(j),(k))]

void fdtd_update_h_cpml(fdtd_grid_t *g, cpml_t *pml) {
    if (!g || !pml) return;
    double cb = g->dt / (MU0 * g->dx);
    double dx = g->dx;

    for (size_t i = 0; i <= g->nx; i++) {
        for (size_t j = 0; j < g->ny; j++) {
            for (size_t k = 0; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                double dEz_dy = EZ(i,j+1,k) - EZ(i,j,k);
                double dEy_dz = EY(i,j,k+1) - EY(i,j,k);
                double termA = cpml_term(dEz_dy, pml->kappa_inv_h_y[j], pml->b_h_y[j], pml->a_h_y[j], dx, &pml->psi_hx_y[idx]);
                double termB = cpml_term(dEy_dz, pml->kappa_inv_h_z[k], pml->b_h_z[k], pml->a_h_z[k], dx, &pml->psi_hx_z[idx]);
                HX(i,j,k) -= cb * (termA - termB);
            }
        }
    }

    for (size_t i = 0; i < g->nx; i++) {
        for (size_t j = 0; j <= g->ny; j++) {
            for (size_t k = 0; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                double dEx_dz = EX(i,j,k+1) - EX(i,j,k);
                double dEz_dx = EZ(i+1,j,k) - EZ(i,j,k);
                double termA = cpml_term(dEx_dz, pml->kappa_inv_h_z[k], pml->b_h_z[k], pml->a_h_z[k], dx, &pml->psi_hy_z[idx]);
                double termB = cpml_term(dEz_dx, pml->kappa_inv_h_x[i], pml->b_h_x[i], pml->a_h_x[i], dx, &pml->psi_hy_x[idx]);
                HY(i,j,k) -= cb * (termA - termB);
            }
        }
    }

    for (size_t i = 0; i < g->nx; i++) {
        for (size_t j = 0; j < g->ny; j++) {
            for (size_t k = 0; k <= g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                double dEy_dx = EY(i+1,j,k) - EY(i,j,k);
                double dEx_dy = EX(i,j+1,k) - EX(i,j,k);
                double termA = cpml_term(dEy_dx, pml->kappa_inv_h_x[i], pml->b_h_x[i], pml->a_h_x[i], dx, &pml->psi_hz_x[idx]);
                double termB = cpml_term(dEx_dy, pml->kappa_inv_h_y[j], pml->b_h_y[j], pml->a_h_y[j], dx, &pml->psi_hz_y[idx]);
                HZ(i,j,k) -= cb * (termA - termB);
            }
        }
    }
}

void fdtd_update_e_cpml(fdtd_grid_t *g, cpml_t *pml) {
    if (!g || !pml) return;
    double ca = g->dt / (EPS0 * g->dx);
    double dx = g->dx;

    for (size_t i = 0; i < g->nx; i++) {
        for (size_t j = 1; j < g->ny; j++) {
            for (size_t k = 1; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                double dHz_dy = HZ(i,j,k) - HZ(i,j-1,k);
                double dHy_dz = HY(i,j,k) - HY(i,j,k-1);
                double termA = cpml_term(dHz_dy, pml->kappa_inv_m_y[j], pml->b_m_y[j], pml->a_m_y[j], dx, &pml->psi_ex_y[idx]);
                double termB = cpml_term(dHy_dz, pml->kappa_inv_m_z[k], pml->b_m_z[k], pml->a_m_z[k], dx, &pml->psi_ex_z[idx]);
                EX(i,j,k) += ca * (termA - termB);
            }
        }
    }

    for (size_t i = 1; i < g->nx; i++) {
        for (size_t j = 0; j < g->ny; j++) {
            for (size_t k = 1; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                double dHx_dz = HX(i,j,k) - HX(i,j,k-1);
                double dHz_dx = HZ(i,j,k) - HZ(i-1,j,k);
                double termA = cpml_term(dHx_dz, pml->kappa_inv_m_z[k], pml->b_m_z[k], pml->a_m_z[k], dx, &pml->psi_ey_z[idx]);
                double termB = cpml_term(dHz_dx, pml->kappa_inv_m_x[i], pml->b_m_x[i], pml->a_m_x[i], dx, &pml->psi_ey_x[idx]);
                EY(i,j,k) += ca * (termA - termB);
            }
        }
    }

    for (size_t i = 1; i < g->nx; i++) {
        for (size_t j = 1; j < g->ny; j++) {
            for (size_t k = 0; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                double dHy_dx = HY(i,j,k) - HY(i-1,j,k);
                double dHx_dy = HX(i,j,k) - HX(i,j-1,k);
                double termA = cpml_term(dHy_dx, pml->kappa_inv_m_x[i], pml->b_m_x[i], pml->a_m_x[i], dx, &pml->psi_ez_x[idx]);
                double termB = cpml_term(dHx_dy, pml->kappa_inv_m_y[j], pml->b_m_y[j], pml->a_m_y[j], dx, &pml->psi_ez_y[idx]);
                EZ(i,j,k) += ca * (termA - termB);
            }
        }
    }
}

#undef EX
#undef EY
#undef EZ
#undef HX
#undef HY
#undef HZ
