#ifndef ENGINE_FDTD_CPML_H
#define ENGINE_FDTD_CPML_H

#include <stddef.h>
#include "grid.h"

struct cpml {
    size_t npml;
    double dx;

    double *kappa_inv_m_x, *b_m_x, *a_m_x;
    double *kappa_inv_m_y, *b_m_y, *a_m_y;
    double *kappa_inv_m_z, *b_m_z, *a_m_z;

    double *kappa_inv_h_x, *b_h_x, *a_h_x;
    double *kappa_inv_h_y, *b_h_y, *a_h_y;
    double *kappa_inv_h_z, *b_h_z, *a_h_z;

    double *psi_hx_y, *psi_hx_z;
    double *psi_hy_z, *psi_hy_x;
    double *psi_hz_x, *psi_hz_y;
    double *psi_ex_y, *psi_ex_z;
    double *psi_ey_z, *psi_ey_x;
    double *psi_ez_x, *psi_ez_y;
};

typedef struct cpml cpml_t;

typedef struct {
    size_t npml;
    double m;
    double ma;
    double kappa_max;
    double alpha_max;
    double r0;
} cpml_params_t;

cpml_params_t cpml_default_params(size_t npml);

cpml_t *cpml_create(const fdtd_grid_t *g, cpml_params_t params);
void cpml_destroy(cpml_t *pml);

void fdtd_update_h_cpml(fdtd_grid_t *g, cpml_t *pml);
void fdtd_update_e_cpml(fdtd_grid_t *g, cpml_t *pml);

static inline double cpml_term(double diff, double kappa_inv, double b, double a,
                                double dx, double *psi)
{
    *psi = b * (*psi) + a * diff;
    return diff * kappa_inv + dx * (*psi);
}

#endif
