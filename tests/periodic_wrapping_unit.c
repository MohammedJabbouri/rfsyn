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

static size_t raw_index(const fdtd_grid_t *g,
                        size_t i, size_t j, size_t k) {
    return (i * (g->ny + 1) + j) * (g->nz + 1) + k;
}

static size_t previous(size_t q, size_t extent) {
    return q == 0 ? extent - 1 : q - 1;
}

static double pattern(size_t i, size_t j, size_t k, unsigned component) {
    size_t value = (31 * i + 17 * j + 7 * k + 5 * component) % 29;
    return ((double)value - 14.0) * 1.0e-4;
}

static void initialize(fdtd_grid_t *g) {
    double *fields[6] = {
        g->ex, g->ey, g->ez, g->hx, g->hy, g->hz
    };

    for (size_t i = 0; i <= g->nx; i++) {
        for (size_t j = 0; j <= g->ny; j++) {
            for (size_t k = 0; k <= g->nz; k++) {
                size_t idx = raw_index(g, i, j, k);

                for (unsigned c = 0; c < 6; c++) {
                    fields[c][idx] = pattern(i, j, k, c);
                }
            }
        }
    }
}

static int active(const fdtd_grid_t *g, int component,
                  size_t i, size_t j, size_t k) {
    if (i >= g->nx || j >= g->ny || k >= g->nz) return 0;

    if (component != 0 && i == 0 && !(g->periodic_axes & 1u)) {
        return 0;
    }

    if (component != 1 && j == 0 && !(g->periodic_axes & 2u)) {
        return 0;
    }

    if (component != 2 && k == 0 && !(g->periodic_axes & 4u)) {
        return 0;
    }

    return 1;
}

// The reference never calls fdtd_index() or the solver's wrapping helper.
static void differences(const fdtd_grid_t *g, int component,
                        size_t i, size_t j, size_t k,
                        double *a, double *b) {
    size_t idx = raw_index(g, i, j, k);

    if (component == 0) {
        size_t jm = previous(j, g->ny);
        size_t km = previous(k, g->nz);

        *a = g->hz[idx] - g->hz[raw_index(g, i, jm, k)];
        *b = g->hy[idx] - g->hy[raw_index(g, i, j, km)];
    } else if (component == 1) {
        size_t im = previous(i, g->nx);
        size_t km = previous(k, g->nz);

        *a = g->hx[idx] - g->hx[raw_index(g, i, j, km)];
        *b = g->hz[idx] - g->hz[raw_index(g, im, j, k)];
    } else {
        size_t im = previous(i, g->nx);
        size_t jm = previous(j, g->ny);

        *a = g->hy[idx] - g->hy[raw_index(g, im, j, k)];
        *b = g->hx[idx] - g->hx[raw_index(g, i, jm, k)];
    }
}

// A fresh CPML object has zero memory, so its first update has this form.
static double stretched(double difference, double inverse_kappa,
                        double a, double dx) {
    return difference * inverse_kappa + dx * (a * difference);
}

static void apply_reference_cpml(const fdtd_grid_t *g, const cpml_t *pml,
                                 int component,
                                 size_t i, size_t j, size_t k,
                                 double *a, double *b) {
    if (component == 0) {
        *a = stretched(*a, pml->kappa_inv_m_y[j], pml->a_m_y[j], g->dx);
        *b = stretched(*b, pml->kappa_inv_m_z[k], pml->a_m_z[k], g->dx);
    } else if (component == 1) {
        *a = stretched(*a, pml->kappa_inv_m_z[k], pml->a_m_z[k], g->dx);
        *b = stretched(*b, pml->kappa_inv_m_x[i], pml->a_m_x[i], g->dx);
    } else {
        *a = stretched(*a, pml->kappa_inv_m_x[i], pml->a_m_x[i], g->dx);
        *b = stretched(*b, pml->kappa_inv_m_y[j], pml->a_m_y[j], g->dx);
    }
}

static int run_case(unsigned axes, int mode) {
    static const field_component_t components[3] = {
        FIELD_EX, FIELD_EY, FIELD_EZ
    };
    static const char *names[4] = {
        "vacuum", "materials", "cpml", "cpml_materials"
    };

    fdtd_grid_t *g = NULL;
    materials_t *mat = NULL;
    cpml_t *pml = NULL;
    double *expected = NULL;
    int result = 1;
    int failures = 0;

    g = fdtd_grid_create(5, 6, 7, 1.0e-3, 0.5);

    if (!g || fdtd_grid_set_periodic(g, axes) != 0) {
        goto cleanup;
    }

    initialize(g);

    size_t total = (g->nx + 1) * (g->ny + 1) * (g->nz + 1);
    double *electric[3] = {g->ex, g->ey, g->ez};

    expected = malloc(3 * total * sizeof(*expected));
    if (!expected) goto cleanup;

    if (mode == 1 || mode == 3) {
        mat = materials_create(g);

        if (!mat ||
            materials_set_box(
                mat, g, 0, g->nx, 0, g->ny, 0, g->nz,
                4.0, 0.02
            ) != 0 ||
            materials_set_component(
                mat, g, FIELD_EY, 1, 1, 1, 9.0, 0.05
            ) != 0) {
            goto cleanup;
        }
    }

    if (mode == 2 || mode == 3) {
        pml = cpml_create(g, cpml_default_params(2));
        if (!pml) goto cleanup;
    }

    for (int c = 0; c < 3; c++) {
        for (size_t i = 0; i <= g->nx; i++) {
            for (size_t j = 0; j <= g->ny; j++) {
                for (size_t k = 0; k <= g->nz; k++) {
                    size_t idx = raw_index(g, i, j, k);
                    size_t slot = (size_t)c * total + idx;
                    double initial = electric[c][idx];

                    expected[slot] = initial;

                    if (!active(g, c, i, j, k)) continue;

                    double a;
                    double b;
                    material_coeff_t coeff = {
                        .Ca = 1.0,
                        .Cb = g->dt / (EPS0 * g->dx)
                    };

                    differences(g, c, i, j, k, &a, &b);

                    if (pml) {
                        apply_reference_cpml(g, pml, c, i, j, k, &a, &b);
                    }

                    if (mat &&
                        materials_get_coefficients(
                            mat, g, components[c], i, j, k, &coeff
                        ) != 0) {
                        goto cleanup;
                    }

                    expected[slot] =
                        coeff.Ca * initial + coeff.Cb * (a - b);
                }
            }
        }
    }

    switch (mode) {
        case 0:
            fdtd_update_e(g);
            break;
        case 1:
            fdtd_update_e_materials(g, mat);
            break;
        case 2:
            fdtd_update_e_cpml(g, pml);
            break;
        case 3:
            fdtd_update_e_cpml_materials(g, pml, mat);
            break;
    }

    for (int c = 0; c < 3; c++) {
        for (size_t i = 0; i <= g->nx; i++) {
            for (size_t j = 0; j <= g->ny; j++) {
                for (size_t k = 0; k <= g->nz; k++) {
                    size_t idx = raw_index(g, i, j, k);
                    double got = electric[c][idx];
                    double want = expected[(size_t)c * total + idx];

                    int matches = active(g, c, i, j, k)
                        ? isfinite(got) &&
                          fabs(got - want) <= 1.0e-12 * (1.0 + fabs(want))
                        : got == want;

                    if (!matches) {
                        if (failures < 8) {
                            fprintf(
                                stderr,
                                "FAIL mask=%u mode=%s component=%d "
                                "(%zu,%zu,%zu): got %.17g expected %.17g\n",
                                axes, names[mode], c,
                                i, j, k, got, want
                            );
                        }

                        failures++;
                    }
                }
            }
        }
    }

    if (failures) {
        fprintf(
            stderr,
            "mask=%u mode=%s: %d mismatches\n",
            axes, names[mode], failures
        );
    }

    result = failures != 0;

cleanup:
    if (result && failures == 0) {
        fprintf(
            stderr,
            "Setup failed: mask=%u mode=%s\n",
            axes, names[mode]
        );
    }

    free(expected);
    cpml_destroy(pml);
    materials_destroy(mat);
    fdtd_grid_destroy(g);
    return result;
}

int main(void) {
    int failed_cases = 0;

    for (unsigned axes = 0; axes < 8; axes++) {
        for (int mode = 0; mode < 4; mode++) {
            failed_cases += run_case(axes, mode);
        }
    }

    if (failed_cases) {
        fprintf(
            stderr,
            "periodic_wrapping_unit: %d of 32 cases failed\n",
            failed_cases
        );
        return EXIT_FAILURE;
    }

    puts("periodic_wrapping_unit: PASS (32 cases, 5 x 6 x 7)");
    return EXIT_SUCCESS;
}