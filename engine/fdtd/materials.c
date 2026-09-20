#include "materials.h"
#include "../core/constants.h"
#include <stdlib.h>

struct materials {
    double *Ca;
    double *Cb;
    size_t total;
};

materials_t *materials_create(const fdtd_grid_t *g) {
    if (!g) return NULL;

    materials_t *mat = malloc(sizeof(materials_t));
    if (!mat) return NULL;

    size_t total = (g->nx + 1) * (g->ny + 1) * (g->nz + 1);
    mat->Ca = malloc(total * sizeof(double));
    mat->Cb = malloc(total * sizeof(double));
    if (!mat->Ca || !mat->Cb) {
        free(mat->Ca);
        free(mat->Cb);
        free(mat);
        return NULL;
    }
    mat->total = total;

    double cb_vacuum = g->dt / (EPS0 * g->dx);
    for (size_t idx = 0; idx < total; idx++) {
        mat->Ca[idx] = 1.0;
        mat->Cb[idx] = cb_vacuum;
    }

    return mat;
}

void materials_destroy(materials_t *mat) {
    if (!mat) return;
    free(mat->Ca);
    free(mat->Cb);
    free(mat);
}

int materials_set_box(materials_t *mat, const fdtd_grid_t *g, size_t i0, size_t i1, size_t j0, size_t j1, size_t k0, size_t k1, double eps_r, double sigma) {

    if (!mat || !g) return -1;
    if (eps_r < 1.0 || sigma < 0.0) return -1;
    if (i1 > g->nx + 1 || j1 > g->ny + 1 || k1 > g->nz + 1) return -1;
    if (i0 >= i1 || j0 >= j1 || k0 >= k1) return -1;

    double eps = eps_r * EPS0;
    double half_sigma_dt_over_eps = (sigma * g->dt) / (2.0 * eps);

    double Ca_val = (1.0 - half_sigma_dt_over_eps) / (1.0 + half_sigma_dt_over_eps);
    double Cb_val = (g->dt / (eps * g->dx)) / (1.0 + half_sigma_dt_over_eps);

    for (size_t i = i0; i < i1; i++) {
        for (size_t j = j0; j < j1; j++) {
            for (size_t k = k0; k < k1; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                mat->Ca[idx] = Ca_val;
                mat->Cb[idx] = Cb_val;
            }
        }
    }

    return 0;
}

#define EX(i,j,k) g->ex[fdtd_index(g,(i),(j),(k))]
#define EY(i,j,k) g->ey[fdtd_index(g,(i),(j),(k))]
#define EZ(i,j,k) g->ez[fdtd_index(g,(i),(j),(k))]
#define HX(i,j,k) g->hx[fdtd_index(g,(i),(j),(k))]
#define HY(i,j,k) g->hy[fdtd_index(g,(i),(j),(k))]
#define HZ(i,j,k) g->hz[fdtd_index(g,(i),(j),(k))]

void fdtd_update_e_materials(fdtd_grid_t *g, const materials_t *mat) {
    if (!g || !mat) return;

    for (size_t i = 0; i < g->nx; i++) {
        for (size_t j = 1; j < g->ny; j++) {
            for (size_t k = 1; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                double curl = (HZ(i,j,k) - HZ(i,j-1,k)) - (HY(i,j,k) - HY(i,j,k-1));
                EX(i,j,k) = mat->Ca[idx] * EX(i,j,k) + mat->Cb[idx] * curl;
            }
        }
    }

    for (size_t i = 1; i < g->nx; i++) {
        for (size_t j = 0; j < g->ny; j++) {
            for (size_t k = 1; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                double curl = (HX(i,j,k) - HX(i,j,k-1)) - (HZ(i,j,k) - HZ(i-1,j,k));
                EY(i,j,k) = mat->Ca[idx] * EY(i,j,k) + mat->Cb[idx] * curl;
            }
        }
    }

    for (size_t i = 1; i < g->nx; i++) {
        for (size_t j = 1; j < g->ny; j++) {
            for (size_t k = 0; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                double curl = (HY(i,j,k) - HY(i-1,j,k)) - (HX(i,j,k) - HX(i,j-1,k));
                EZ(i,j,k) = mat->Ca[idx] * EZ(i,j,k) + mat->Cb[idx] * curl;
            }
        }
    }
}

void fdtd_update_e_cpml_materials(fdtd_grid_t *g, cpml_t *pml, const materials_t *mat) {
    if (!g || !pml || !mat) return;
    double dx = g->dx;

    for (size_t i = 0; i < g->nx; i++) {
        for (size_t j = 1; j < g->ny; j++) {
            for (size_t k = 1; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                double dHz_dy = HZ(i,j,k) - HZ(i,j-1,k);
                double dHy_dz = HY(i,j,k) - HY(i,j,k-1);
                double termA = cpml_term(dHz_dy, pml->kappa_inv_m_y[j], pml->b_m_y[j], pml->a_m_y[j], dx, &pml->psi_ex_y[idx]);
                double termB = cpml_term(dHy_dz, pml->kappa_inv_m_z[k], pml->b_m_z[k], pml->a_m_z[k], dx, &pml->psi_ex_z[idx]);
                EX(i,j,k) = mat->Ca[idx] * EX(i,j,k) + mat->Cb[idx] * (termA - termB);
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
                EY(i,j,k) = mat->Ca[idx] * EY(i,j,k) + mat->Cb[idx] * (termA - termB);
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
                EZ(i,j,k) = mat->Ca[idx] * EZ(i,j,k) + mat->Cb[idx] * (termA - termB);
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