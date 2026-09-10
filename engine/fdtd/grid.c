#include "grid.h"
#include "../core/constants.h"
#include <stdlib.h>
#include <math.h>

fdtd_grid_t *fdtd_grid_create(size_t nx, size_t ny, size_t nz, double dx, double courant_safety_factor) {
    if (dx <= 0.0 || courant_safety_factor <= 0.0) return NULL;
    if (nx > SIZE_MAX - 1 || ny > SIZE_MAX - 1 || nz > SIZE_MAX - 1) return NULL;
    
    size_t nx1 = nx + 1;
    size_t ny1 = ny + 1;
    size_t nz1 = nz + 1;
    
    if (nx1 > SIZE_MAX / ny1 || (nx1 * ny1) > SIZE_MAX / nz1) return NULL;

    fdtd_grid_t *g = malloc(sizeof(fdtd_grid_t));
    if (!g) return NULL;
    
    g->nx = nx;
    g->ny = ny;
    g->nz = nz;
    g->dx = dx;

    double dt_max = dx / (C0 * sqrt(3.0));
    g->dt = courant_safety_factor * dt_max;
    
    size_t total = nx1 * ny1 * nz1;
    
    g->ex = calloc(total, sizeof(double));
    g->ey = calloc(total, sizeof(double));
    g->ez = calloc(total, sizeof(double));
    g->hx = calloc(total, sizeof(double));
    g->hy = calloc(total, sizeof(double));
    g->hz = calloc(total, sizeof(double));

    if (!g->ex || !g->ey || !g->ez || !g->hx || !g->hy || !g->hz) {
        fdtd_grid_destroy(g);
        return NULL;
    }

    return g;
}

void fdtd_grid_destroy(fdtd_grid_t *g) {
    if (!g) return;
    free(g->ex); free(g->ey); free(g->ez);
    free(g->hx); free(g->hy); free(g->hz);
    free(g);
}

#define EX(i,j,k) g->ex[fdtd_index(g,(i),(j),(k))]
#define EY(i,j,k) g->ey[fdtd_index(g,(i),(j),(k))]
#define EZ(i,j,k) g->ez[fdtd_index(g,(i),(j),(k))]
#define HX(i,j,k) g->hx[fdtd_index(g,(i),(j),(k))]
#define HY(i,j,k) g->hy[fdtd_index(g,(i),(j),(k))]
#define HZ(i,j,k) g->hz[fdtd_index(g,(i),(j),(k))]

void fdtd_update_h(fdtd_grid_t *g) {
    double cb = g->dt / (MU0 * g->dx);

    for (size_t i = 0; i <= g->nx; i++) {
        for (size_t j = 0; j < g->ny; j++) {
            for (size_t k = 0; k < g->nz; k++) {
                HX(i,j,k) -= cb * ((EZ(i,j+1,k) - EZ(i,j,k)) - (EY(i,j,k+1) - EY(i,j,k)));
            }
        }
    }

    for (size_t i = 0; i < g->nx; i++) {
        for (size_t j = 0; j <= g->ny; j++) {
            for (size_t k = 0; k < g->nz; k++) {
                HY(i,j,k) -= cb * ((EX(i,j,k+1) - EX(i,j,k)) - (EZ(i+1,j,k) - EZ(i,j,k)));
            }
        }
    }

    for (size_t i = 0; i < g->nx; i++) {
        for (size_t j = 0; j < g->ny; j++) {
            for (size_t k = 0; k <= g->nz; k++) {
                HZ(i,j,k) -= cb * ((EY(i+1,j,k) - EY(i,j,k)) - (EX(i,j+1,k) - EX(i,j,k)));
            }
        }
    }
}

void fdtd_update_e(fdtd_grid_t *g) {
    double ca = g->dt / (EPS0 * g->dx);

    for (size_t i = 0; i < g->nx; i++) {
        for (size_t j = 1; j < g->ny; j++) {
            for (size_t k = 1; k < g->nz; k++) {
                EX(i,j,k) += ca * ((HZ(i,j,k) - HZ(i,j-1,k)) - (HY(i,j,k) - HY(i,j,k-1)));
            }
        }
    }

    for (size_t i = 1; i < g->nx; i++) {
        for (size_t j = 0; j < g->ny; j++) {
            for (size_t k = 1; k < g->nz; k++) {
                EY(i,j,k) += ca * ((HX(i,j,k) - HX(i,j,k-1)) - (HZ(i,j,k) - HZ(i-1,j,k)));
            }
        }
    }

    for (size_t i = 1; i < g->nx; i++) {
        for (size_t j = 1; j < g->ny; j++) {
            for (size_t k = 0; k < g->nz; k++) {
                EZ(i,j,k) += ca * ((HY(i,j,k) - HY(i-1,j,k)) - (HX(i,j,k) - HX(i,j-1,k)));
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
