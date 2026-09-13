#ifndef ENGINE_FDTD_CPML_H
#define ENGINE_FDTD_CPML_H

#include <stddef.h>
#include "grid.h"

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

#endif
