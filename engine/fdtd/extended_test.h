#ifndef RFSYN_EXTENDED_TEST_H
#define RFSYN_EXTENDED_TEST_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <errno.h>

#include "grid.h"
#include "cpml.h"
#include "materials.h"
#include "../core/constants.h"

#ifndef RFSYN_VALIDATION_REFINEMENT
#define RFSYN_VALIDATION_REFINEMENT 1
#endif

#define XT_PI 3.1415926535897932384626433832795

typedef struct {
    fdtd_grid_t *g;
    cpml_t *p;
    materials_t *m;
} xt_sim;

static inline void xt_destroy(xt_sim *s)
{
    if (!s) return;

    materials_destroy(s->m);
    cpml_destroy(s->p);
    fdtd_grid_destroy(s->g);

    memset(s, 0, sizeof(*s));
}

static inline int xt_create(
    xt_sim *s,
    size_t nx,
    size_t ny,
    size_t nz,
    double dx,
    unsigned periodic_axes,
    size_t npml)
{
    memset(s, 0, sizeof(*s));

    s->g = fdtd_grid_create(nx, ny, nz, dx, 0.99);

    if (!s->g) {
        xt_destroy(s);
        return -1;
    }

    if (fdtd_grid_set_periodic(s->g, periodic_axes) != 0) {
        xt_destroy(s);
        return -1;
    }

    if (npml > 0) {
        s->p = cpml_create(s->g, cpml_default_params(npml));

        if (!s->p) {
            xt_destroy(s);
            return -1;
        }
    }

    s->m = materials_create(s->g);

    if (!s->m) {
        xt_destroy(s);
        return -1;
    }

    return 0;
}

static inline void xt_step(xt_sim *s)
{
    if (s->p) {
        fdtd_update_h_cpml(s->g, s->p);
        fdtd_update_e_cpml_materials(s->g, s->p, s->m);
    } else {
        fdtd_update_h(s->g);
        fdtd_update_e_materials(s->g, s->m);
    }
}

static inline int xt_medium(xt_sim *s, double eps_r)
{
    if (!s || !s->g || !s->m) {
        return -1;
    }

    const fdtd_grid_t *g = s->g;

    size_t ix = g->nx + ((g->periodic_axes & 1u) ? 0u : 1u);
    size_t jy = g->ny + ((g->periodic_axes & 2u) ? 0u : 1u);
    size_t kz = g->nz + ((g->periodic_axes & 4u) ? 0u : 1u);

    return materials_set_box(
        s->m, g,
        0, ix,
        0, jy,
        0, kz,
        eps_r, 0.0
    );
}

static inline int xt_interface(
    xt_sim *s,
    size_t i_interface,
    double eps_left,
    double eps_right)
{
    if (!s || !s->g || !s->m ||
        i_interface == 0 ||
        i_interface >= s->g->nx) {
        return -1;
    }

    if (xt_medium(s, eps_left) != 0) {
        return -1;
    }

    const fdtd_grid_t *g = s->g;

    size_t ix = g->nx + ((g->periodic_axes & 1u) ? 0u : 1u);
    size_t jy = g->ny + ((g->periodic_axes & 2u) ? 0u : 1u);
    size_t kz = g->nz + ((g->periodic_axes & 4u) ? 0u : 1u);

    if (materials_set_box(
            s->m, g,
            i_interface, ix,
            0, jy,
            0, kz,
            eps_right, 0.0) != 0) {
        return -1;
    }

    double eps_tangent = 0.5 * (eps_left + eps_right);
    size_t planes[2] = {i_interface, 0};
    size_t plane_count = (g->periodic_axes & 1u) ? 2u : 1u;

    for (size_t p = 0; p < plane_count; p++) {
        for (size_t j = 0; j < g->ny; j++) {
            for (size_t k = 0; k < g->nz; k++) {
                if (materials_set_component(
                        s->m, g, FIELD_EY,
                        planes[p], j, k,
                        eps_tangent, 0.0) != 0) {
                    return -1;
                }

                if (materials_set_component(
                        s->m, g, FIELD_EZ,
                        planes[p], j, k,
                        eps_tangent, 0.0) != 0) {
                    return -1;
                }
            }
        }
    }

    return 0;
}

static inline double *xt_e(fdtd_grid_t *g, int tm)
{
    return tm ? g->ey : g->ez;
}

static inline int xt_finite_complex(double complex z)
{
    return isfinite(creal(z)) && isfinite(cimag(z));
}

static inline double xt_db(double amplitude)
{
    return amplitude > 0.0 ? 20.0 * log10(amplitude) : -INFINITY;
}

static inline double xt_kx(
    double omega_dt,
    double courant,
    double eps_r,
    double ky)
{
    double z =
        eps_r * pow(sin(omega_dt / 2.0) / courant, 2.0)
        - pow(sin(ky / 2.0), 2.0);

    if (!(z > 0.0 && z < 1.0)) return NAN;

    return 2.0 * asin(sqrt(z));
}

static inline void xt_modes(
    double complex e1,
    double complex e2,
    double phase,
    double complex *forward,
    double complex *backward)
{
    double complex positive = cexp(I * phase);
    double complex negative = conj(positive);
    double complex denominator = positive - negative;

    *forward = (e2 - e1 * negative) / denominator;
    *backward = (e1 * positive - e2) / denominator;
}

static inline int xt_parse_full(int argc, char **argv, int *full)
{
    *full = 0;

    for (int a = 1; a < argc; a++) {
        if (strcmp(argv[a], "--full") == 0) {
            *full = 1;
        } else {
            return -1;
        }
    }

    return 0;
}

static inline int xt_parse_stability(
    int argc,
    char **argv,
    size_t *steps,
    int *full)
{
    *steps = 100000;
    *full = 0;

    for (int a = 1; a < argc; a++) {
        if (strcmp(argv[a], "--full") == 0) {
            *full = 1;
        } else if (strcmp(argv[a], "--steps") == 0 && a + 1 < argc) {
            char *end = NULL;
            errno = 0;

            unsigned long long value = strtoull(argv[++a], &end, 10);

            if (errno != 0 ||
                end == argv[a] ||
                *end != '\0' ||
                value < 10000 ||
                value > 10000000) {
                return -1;
            }

            *steps = (size_t)value;
        } else {
            return -1;
        }
    }

    return 0;
}

#endif