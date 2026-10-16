#ifndef ENGINE_FDTD_GRID_H
#define ENGINE_FDTD_GRID_H

#include <stddef.h>

typedef struct {
    size_t nx, ny, nz;
    double dx;
    double dt;
    unsigned periodic_axes;
    double *ex, *ey, *ez;
    double *hx, *hy, *hz;
} fdtd_grid_t;

fdtd_grid_t *fdtd_grid_create(size_t nx, size_t ny, size_t nz, double dx, double courant_safety_factor);

// Set before creating materials/CPML or initializing fields. Bits: x=1,y=2,z=4.

int fdtd_grid_set_periodic(fdtd_grid_t *g, unsigned axes);

void fdtd_grid_destroy(fdtd_grid_t *g);

static inline size_t fdtd_index(const fdtd_grid_t *g, size_t i, size_t j, size_t k) {
    if (g->periodic_axes & 1u) i %= g->nx;
    if (g->periodic_axes & 2u) j %= g->ny;
    if (g->periodic_axes & 4u) k %= g->nz;
    return (i * (g->ny + 1) + j) * (g->nz + 1) + k;
}

void fdtd_update_h(fdtd_grid_t *g);
void fdtd_update_e(fdtd_grid_t *g);

#endif