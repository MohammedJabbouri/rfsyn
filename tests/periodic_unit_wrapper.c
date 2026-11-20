//
// Created by root on 11/20/26.
//

#include "../engine/fdtd/grid.h"
#include "../engine/fdtd/materials.h"
#include "../engine/fdtd/cpml.h"
#include "../engine/core/constants.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Regression test: prevents size_t underflow (SIZE_MAX % n) during periodic backward-neighbor wrapping on a 5x6x7 grid.

#define NX 5
#define NY 6
#define NZ 7
#define CELL 1.0e-3
#define COURANT 0.5

static int failures = 0;

static void report(const char *name, int comp, size_t i, size_t j, size_t k,
                   double got, double want) {
    printf("FAIL %s c%d(%zu,%zu,%zu): got %.17g want %.17g\n",
           name, comp, i, j, k, got, want);
    failures++;
}

static size_t raw_index(const fdtd_grid_t *g, size_t i, size_t j, size_t k) {
    return (i * (g->ny + 1) + j) * (g->nz + 1) + k;
}

static size_t prev_coord(size_t q, size_t n, unsigned periodic) {
    return periodic ? (q == 0 ? n - 1 : q - 1) : q - 1;
}

static size_t lo_x(unsigned axes) { return (axes & 1u) ? 0u : 1u; }
static size_t lo_y(unsigned axes) { return (axes & 2u) ? 0u : 1u; }
static size_t lo_z(unsigned axes) { return (axes & 4u) ? 0u : 1u; }

static void expect_close_probe(const char *name, const fdtd_grid_t *g,
                               int comp, size_t i, size_t j, size_t k,
                               double got, double want) {
    (void)g;
    if (fabs(got - want) > 1.0e-12 * (1.0 + fabs(want)))
        report(name, comp, i, j, k, got, want);
}

static int close_enough(double got, double want) {
    return fabs(got - want) <= 1.0e-12 * (1.0 + fabs(want));
}

static double pattern(size_t i, size_t j, size_t k, int field_id) {
    unsigned m = (unsigned)((i * 31u + j * 17u + k * 7u
                             + (size_t)field_id * 5u) % 29u);
    return ((double)m - 14.0) * 1.0e-4;
}

static void fill_fields(fdtd_grid_t *g) {
    for (size_t i = 0; i <= g->nx; i++) {
        for (size_t j = 0; j <= g->ny; j++) {
            for (size_t k = 0; k <= g->nz; k++) {
                size_t idx = raw_index(g, i, j, k);
                g->ex[idx] = pattern(i, j, k, 0);
                g->ey[idx] = pattern(i, j, k, 1);
                g->ez[idx] = pattern(i, j, k, 2);
                g->hx[idx] = pattern(i, j, k, 3);
                g->hy[idx] = pattern(i, j, k, 4);
                g->hz[idx] = pattern(i, j, k, 5);
            }
        }
    }
}

static size_t total_cells(const fdtd_grid_t *g) {
    return (g->nx + 1) * (g->ny + 1) * (g->nz + 1);
}

static double *copy_e_fields(const fdtd_grid_t *g) {
    size_t total = total_cells(g);
    double *snap = malloc(3 * total * sizeof(double));
    if (!snap) { fprintf(stderr, "allocation failed\n"); exit(2); }
    memcpy(snap, g->ex, total * sizeof(double));
    memcpy(snap + total, g->ey, total * sizeof(double));
    memcpy(snap + 2 * total, g->ez, total * sizeof(double));
    return snap;
}

static void curl_terms(const fdtd_grid_t *g, unsigned axes, int comp, size_t i, size_t j, size_t k, double *term_a, double *term_b) {
    size_t pi = prev_coord(i, g->nx, axes & 1u);
    size_t pj = prev_coord(j, g->ny, axes & 2u);
    size_t pk = prev_coord(k, g->nz, axes & 4u);

    if (comp == 0) {
        *term_a = g->hz[raw_index(g, i, j, k)] - g->hz[raw_index(g, i, pj, k)];
        *term_b = g->hy[raw_index(g, i, j, k)] - g->hy[raw_index(g, i, j, pk)];
    } else if (comp == 1) {
        *term_a = g->hx[raw_index(g, i, j, k)] - g->hx[raw_index(g, i, j, pk)];
        *term_b = g->hz[raw_index(g, i, j, k)] - g->hz[raw_index(g, pi, j, k)];
    } else {
        *term_a = g->hy[raw_index(g, i, j, k)] - g->hy[raw_index(g, pi, j, k)];
        *term_b = g->hx[raw_index(g, i, j, k)] - g->hx[raw_index(g, i, pj, k)];
    }
}

static void domain_bounds(const fdtd_grid_t *g, unsigned axes, int comp, size_t *i0, size_t *j0, size_t *k0) {
    (void)g;
    if (comp == 0)      { *i0 = 0;          *j0 = lo_y(axes); *k0 = lo_z(axes); }
    else if (comp == 1) { *i0 = lo_x(axes); *j0 = 0;          *k0 = lo_z(axes); }
    else                { *i0 = lo_x(axes); *j0 = lo_y(axes); *k0 = 0;          }
}

static int updated_position(const fdtd_grid_t *g, unsigned axes, int comp,
                            size_t i, size_t j, size_t k) {
    size_t i0, j0, k0;
    domain_bounds(g, axes, comp, &i0, &j0, &k0);
    return i0 <= i && i < g->nx && j0 <= j && j < g->ny && k0 <= k && k < g->nz;
}

// fdtd_update_e()
static void reference_vacuum(const fdtd_grid_t *g, unsigned axes, const double *snap, double *ref, size_t total) {
    double ca = g->dt / (EPS0 * g->dx);
    memcpy(ref, snap, 3 * total * sizeof(double));

    for (int c = 0; c < 3; c++) {
        size_t i0, j0, k0;
        domain_bounds(g, axes, c, &i0, &j0, &k0);
        for (size_t i = i0; i < g->nx; i++) {
            for (size_t j = j0; j < g->ny; j++) {
                for (size_t k = k0; k < g->nz; k++) {
                    double ta, tb;
                    size_t idx = raw_index(g, i, j, k);
                    curl_terms(g, axes, c, i, j, k, &ta, &tb);
                    ref[(size_t)c * total + idx] =
                        snap[(size_t)c * total + idx] + ca * (ta - tb);
                }
            }
        }
    }
}

// fdtd_update_e_materials()
static void reference_materials(const fdtd_grid_t *g, const materials_t *mat, unsigned axes, const double *snap, double *ref, size_t total) {
    static const field_component_t comps[3] = {FIELD_EX, FIELD_EY, FIELD_EZ};
    memcpy(ref, snap, 3 * total * sizeof(double));

    for (int c = 0; c < 3; c++) {
        size_t i0, j0, k0;
        domain_bounds(g, axes, c, &i0, &j0, &k0);
        for (size_t i = i0; i < g->nx; i++) {
            for (size_t j = j0; j < g->ny; j++) {
                for (size_t k = k0; k < g->nz; k++) {
                    double ta, tb;
                    size_t idx = raw_index(g, i, j, k);
                    material_coeff_t coeff;
                    if (materials_get_coefficients(mat, g, comps[c], i, j, k, &coeff) != 0) {
                        fprintf(stderr, "coefficient lookup failed\n");
                        exit(2);
                    }
                    curl_terms(g, axes, c, i, j, k, &ta, &tb);
                    ref[(size_t)c * total + idx] =
                        coeff.Ca * snap[(size_t)c * total + idx]
                        + coeff.Cb * (ta - tb);
                }
            }
        }
    }
}

static void cpml_terms(fdtd_grid_t *g, cpml_t *pml, unsigned axes, int comp,
                       size_t i, size_t j, size_t k,
                       double *term_a, double *term_b) {
    double ta, tb;
    size_t idx = raw_index(g, i, j, k);
    curl_terms(g, axes, comp, i, j, k, &ta, &tb);

    if (comp == 0) {
        *term_a = cpml_term(ta, pml->kappa_inv_m_y[j], pml->b_m_y[j], pml->a_m_y[j], g->dx, &pml->psi_ex_y[idx]);
        *term_b = cpml_term(tb, pml->kappa_inv_m_z[k], pml->b_m_z[k], pml->a_m_z[k], g->dx, &pml->psi_ex_z[idx]);
    } else if (comp == 1) {
        *term_a = cpml_term(ta, pml->kappa_inv_m_z[k], pml->b_m_z[k], pml->a_m_z[k], g->dx, &pml->psi_ey_z[idx]);
        *term_b = cpml_term(tb, pml->kappa_inv_m_x[i], pml->b_m_x[i], pml->a_m_x[i], g->dx, &pml->psi_ey_x[idx]);
    } else {
        *term_a = cpml_term(ta, pml->kappa_inv_m_x[i], pml->b_m_x[i], pml->a_m_x[i], g->dx, &pml->psi_ez_x[idx]);
        *term_b = cpml_term(tb, pml->kappa_inv_m_y[j], pml->b_m_y[j], pml->a_m_y[j], g->dx, &pml->psi_ez_y[idx]);
    }
}

static void reference_cpml(fdtd_grid_t *g, cpml_t *pml_ref, const materials_t *mat, unsigned axes, const double *snap, double *ref, size_t total) {
    static const field_component_t comps[3] = {FIELD_EX, FIELD_EY, FIELD_EZ};
    double ca = g->dt / (EPS0 * g->dx);
    memcpy(ref, snap, 3 * total * sizeof(double));

    for (int c = 0; c < 3; c++) {
        size_t i0, j0, k0;
        domain_bounds(g, axes, c, &i0, &j0, &k0);
        for (size_t i = i0; i < g->nx; i++) {
            for (size_t j = j0; j < g->ny; j++) {
                for (size_t k = k0; k < g->nz; k++) {
                    double ta, tb;
                    size_t idx = raw_index(g, i, j, k);
                    cpml_terms(g, pml_ref, axes, c, i, j, k, &ta, &tb);
                    if (mat) {
                        material_coeff_t coeff;
                        if (materials_get_coefficients(mat, g, comps[c],
                                                       i, j, k, &coeff) != 0) {
                            fprintf(stderr, "coefficient lookup failed\n");
                            exit(2);
                        }
                        ref[(size_t)c * total + idx] =
                            coeff.Ca * snap[(size_t)c * total + idx]
                            + coeff.Cb * (ta - tb);
                    } else {
                        ref[(size_t)c * total + idx] =
                            snap[(size_t)c * total + idx] + ca * (ta - tb);
                    }
                }
            }
        }
    }
}

static void compare_result(const char *name, const fdtd_grid_t *g,
                           unsigned axes, const double *snap,
                           const double *ref, size_t total) {
    const double *actual[3] = {g->ex, g->ey, g->ez};
    for (int c = 0; c < 3; c++) {
        for (size_t i = 0; i <= g->nx; i++) {
            for (size_t j = 0; j <= g->ny; j++) {
                for (size_t k = 0; k <= g->nz; k++) {
                    size_t idx = raw_index(g, i, j, k);
                    size_t slot = (size_t)c * total + idx;
                    if (updated_position(g, axes, c, i, j, k)) {
                        if (!close_enough(actual[c][idx], ref[slot]))
                            report(name, c, i, j, k, actual[c][idx], ref[slot]);
                    } else {
                        if (!close_enough(actual[c][idx], snap[slot]))
                            report(name, c, i, j, k, actual[c][idx], snap[slot]);
                    }
                }
            }
        }
    }
}

static void run_reference_case(unsigned axes, int mode, const char *name) {
    fdtd_grid_t *g = fdtd_grid_create(NX, NY, NZ, CELL, COURANT);
    materials_t *mat = NULL;
    cpml_t *pml = NULL, *pml_ref = NULL;
    size_t total;
    double *snap, *ref;

    if (!g) { fprintf(stderr, "grid creation failed\n"); exit(2); }
    if (fdtd_grid_set_periodic(g, axes) != 0) {
        fprintf(stderr, "set_periodic failed\n"); exit(2);
    }
    fill_fields(g);
    total = total_cells(g);
    snap = copy_e_fields(g);
    ref = copy_e_fields(g);

    if (mode == 1 || mode == 3) {
        material_id_t slabs;
        mat = materials_create(g);
        if (!mat) { fprintf(stderr, "materials creation failed\n"); exit(2); }
        if (materials_register(mat, 4.0, 0.02, &slabs) != 0
            || materials_assign_box(mat, g, 1, 4, 0, 6, 2, 5, slabs) != 0
            || materials_set_component(mat, g, FIELD_EY, 0, 0, 0,
                                       6.0, 0.1) != 0) {
            fprintf(stderr, "material setup failed\n"); exit(2);
        }
    }
    if (mode == 2 || mode == 3) {
        pml = cpml_create(g, cpml_default_params(2));
        pml_ref = cpml_create(g, cpml_default_params(2));
        if (!pml || !pml_ref) {
            fprintf(stderr, "cpml creation failed\n"); exit(2);
        }
    }

    switch (mode) {
        case 0:
            reference_vacuum(g, axes, snap, ref, total);
            fdtd_update_e(g);
            break;
        case 1:
            reference_materials(g, mat, axes, snap, ref, total);
            fdtd_update_e_materials(g, mat);
            break;
        case 2:
            reference_cpml(g, pml_ref, NULL, axes, snap, ref, total);
            fdtd_update_e_cpml(g, pml);
            break;
        default:
            reference_cpml(g, pml_ref, mat, axes, snap, ref, total);
            fdtd_update_e_cpml_materials(g, pml, mat);
            break;
    }

    compare_result(name, g, axes, snap, ref, total);

    cpml_destroy(pml_ref);
    cpml_destroy(pml);
    materials_destroy(mat);
    free(ref);
    free(snap);
    fdtd_grid_destroy(g);
}

static void probe_x_wrap(void) {
    fdtd_grid_t *g = fdtd_grid_create(NX, NY, NZ, CELL, COURANT);
    double ca;
    if (!g || fdtd_grid_set_periodic(g, 1u) != 0) exit(2);
    ca = g->dt / (EPS0 * g->dx);
    g->hy[raw_index(g, NX - 1, 2, 3)] = 1.0e-3;
    fdtd_update_e(g);
    expect_close_probe("x_wrap_seam", g, 2, 0, 2, 3, g->ez[raw_index(g, 0, 2, 3)], -ca * 1.0e-3);
    expect_close_probe("x_wrap_far", g, 2, NX - 1, 2, 3, g->ez[raw_index(g, NX - 1, 2, 3)], ca * 1.0e-3);
    fdtd_grid_destroy(g);
}

static void probe_y_wrap(void) {
    fdtd_grid_t *g = fdtd_grid_create(NX, NY, NZ, CELL, COURANT);
    double ca;
    if (!g || fdtd_grid_set_periodic(g, 2u) != 0) exit(2);
    ca = g->dt / (EPS0 * g->dx);
    g->hz[raw_index(g, 2, NY - 1, 3)] = 1.0e-3;
    fdtd_update_e(g);
    expect_close_probe("y_wrap_seam", g, 0, 2, 0, 3, g->ex[raw_index(g, 2, 0, 3)], -ca * 1.0e-3);
    expect_close_probe("y_wrap_far", g, 0, 2, NY - 1, 3, g->ex[raw_index(g, 2, NY - 1, 3)], ca * 1.0e-3);
    fdtd_grid_destroy(g);
}

static void probe_z_wrap(void) {
    fdtd_grid_t *g = fdtd_grid_create(NX, NY, NZ, CELL, COURANT);
    double ca;
    if (!g || fdtd_grid_set_periodic(g, 4u) != 0) exit(2);
    ca = g->dt / (EPS0 * g->dx);
    g->hy[raw_index(g, 2, 3, NZ - 1)] = 1.0e-3;
    fdtd_update_e(g);
    expect_close_probe("z_wrap_seam", g, 0, 2, 3, 0, g->ex[raw_index(g, 2, 3, 0)], ca * 1.0e-3);
    expect_close_probe("z_wrap_far", g, 0, 2, 3, NZ - 1, g->ex[raw_index(g, 2, 3, NZ - 1)], -ca * 1.0e-3);
    fdtd_grid_destroy(g);
}

int main(void) {
    static const char *names[4] = {
        "vacuum", "materials", "cpml", "cpml_materials"
    };

    probe_x_wrap();
    probe_y_wrap();
    probe_z_wrap();

    for (unsigned axes = 0; axes < 8; axes++) {
        for (int mode = 0; mode < 4; mode++) {
            run_reference_case(axes, mode, names[mode]);
        }
    }

    if (failures) {
        printf("periodic_wrapping_unit: %d check(s) FAILED\n", failures);
        return 1;
    }
    puts("periodic_wrapping_unit: PASS");
    return 0;
}