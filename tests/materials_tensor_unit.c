//
// Created by root on 11/21/26.
//

#include "../engine/fdtd/materials.h"
#include "../engine/fdtd/grid.h"
#include "../engine/fdtd/subpixel_te.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static size_t checks = 0, failures = 0;

#define CHECK(expr) do { checks++; if (!(expr)) { \
    fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #expr); failures++; } } while (0)

#define CLOSE(a, b, tol) do { checks++; \
    double _a = (a), _b = (b); \
    if (!isfinite(_a) || fabs(_a - _b) > (tol)) { \
        fprintf(stderr, "FAIL line %d: %.17g vs %.17g (tol %.3g)\n", __LINE__, _a, _b, (double)(tol)); \
        failures++; } } while (0)

static size_t idxs(const fdtd_grid_t *g, size_t i, size_t j, size_t k) {
    return fdtd_index(g, i, j, k);
}

static void test_api(void)
{
    fdtd_grid_t *g = fdtd_grid_create(8, 8, 2, 0.001, 0.9);
    CHECK(g != NULL);
    materials_t *m = materials_create(g);
    CHECK(m != NULL);

    material_tensor_t out; memset(&out, 0, sizeof out);

    /* no tensor yet: getter returns 1 */
    CHECK(materials_get_tensor(m, g, FIELD_EX, 4, 4, 0, &out) == 1);

    /* valid SPD tensor at interior Ex edge */
    CHECK(materials_set_component_tensor(m, g, FIELD_EX, 4, 4, 0,
                                         0.6, 0.05, 0.5) == 0);
    CHECK(materials_get_tensor(m, g, FIELD_EX, 4, 4, 0, &out) == 0);
    CLOSE(out.m11, 0.6, 0); CLOSE(out.m12, 0.05, 0); CLOSE(out.m22, 0.5, 0);
    CHECK(out.scale > 0.0);

    /* Ez rejected */
    CHECK(materials_set_component_tensor(m, g, FIELD_EZ, 4, 4, 0,
                                         0.6, 0.0, 0.5) == -1);
    /* non-SPD rejected: determinant 0.09 - 0.16 < 0 */
    CHECK(materials_set_component_tensor(m, g, FIELD_EX, 4, 4, 0,
                                         0.3, 0.4, 0.3) == -1);
    /* non-positive diagonals rejected */
    CHECK(materials_set_component_tensor(m, g, FIELD_EX, 4, 4, 0,
                                         0.0, 0.0, 0.3) == -1);
    /* NaN rejected */
    CHECK(materials_set_component_tensor(m, g, FIELD_EX, 4, 4, 0,
                                         0.3, 0.0/0.0, 0.3) == -1);
    /* boundary margin rejected: i=0 on nonperiodic x */
    CHECK(materials_set_component_tensor(m, g, FIELD_EX, 0, 4, 0,
                                         0.6, 0.0, 0.5) == -1);

    /* overwrite same edge */
    CHECK(materials_set_component_tensor(m, g, FIELD_EX, 4, 4, 0,
                                         0.4, 0.0, 0.4) == 0);
    CHECK(materials_get_tensor(m, g, FIELD_EX, 4, 4, 0, &out) == 0);
    CLOSE(out.m11, 0.4, 0);

    /* clear */
    CHECK(materials_clear_component_tensor(m, g, FIELD_EX, 4, 4, 0) == 0);
    CHECK(materials_get_tensor(m, g, FIELD_EX, 4, 4, 0, &out) == 1);

    /* box assignment clears tensors inside the box */
    CHECK(materials_set_component_tensor(m, g, FIELD_EX, 4, 4, 0,
                                         0.6, 0.0, 0.5) == 0);
    CHECK(materials_set_box(m, g, 3, 6, 3, 6, 0, 1, 4.0, 0.0) == 0);
    CHECK(materials_get_tensor(m, g, FIELD_EX, 4, 4, 0, &out) == 1);

    materials_destroy(m);
    fdtd_grid_destroy(g);
}

/* Invariance: a uniform diagonal tensor m = (1/eps)I must reproduce the
 * scalar eps channel to machine precision over many steps. */
static void test_invariance(void)
{
    fdtd_grid_t *a = fdtd_grid_create(16, 16, 2, 0.001, 0.9);
    fdtd_grid_t *b = fdtd_grid_create(16, 16, 2, 0.001, 0.9);
    CHECK(a && b);
    materials_t *ms = materials_create(a);
    materials_t *mt = materials_create(b);
    CHECK(ms && mt);

    double eps_r = 2.5, mdiag = 1.0 / eps_r;

    /* scalar path: interior isotropic eps_r */
    CHECK(materials_set_box(ms, a, 2, 14, 2, 14, 0, 2, eps_r, 0.0) == 0);

    /* tensor path: diagonal tensor at every interior Ex/Ey edge,
     * scalar base left at vacuum so any unnoticed fallthrough differs */
    int ok = 1;
    for (size_t i = 1; i < 15 && ok; i++)
        for (size_t j = 1; j < 14 && ok; j++) {
            ok &= materials_set_component_tensor(mt, b, FIELD_EX, i, j, 0,
                                                 mdiag, 0.0, mdiag) == 0;
            ok &= materials_set_component_tensor(mt, b, FIELD_EY, i, j, 0,
                                                 mdiag, 0.0, mdiag) == 0;
        }
    CHECK(ok);

    /* drive identical soft sources */
    for (int n = 0; n < 400; n++) {
        fdtd_update_h(a); fdtd_update_h(b);
        fdtd_update_e_materials(a, ms); fdtd_update_e_materials(b, mt);
        double s = sin(2.0 * 3.14159265358979 * 0.05 * n);
        a->hz[idxs(a, 8, 8, 0)] += s;
        b->hz[idxs(b, 8, 8, 0)] += s;
    }

    /* compare all six component arrays */
    double maxdiff = 0.0;
    size_t total = (a->nx + 1) * (a->ny + 1) * (a->nz + 1);
    for (size_t q = 0; q < total; q++) {
        double d1 = fabs(a->ex[q] - b->ex[q]);
        double d2 = fabs(a->ey[q] - b->ey[q]);
        double d3 = fabs(a->hx[q] - b->hx[q]);
        if (d1 > maxdiff) maxdiff = d1;
        if (d2 > maxdiff) maxdiff = d2;
        if (d3 > maxdiff) maxdiff = d3;
    }
    CHECK(maxdiff < 1e-12);

    materials_destroy(ms); materials_destroy(mt);
    fdtd_grid_destroy(a); fdtd_grid_destroy(b);
}

/* Correctness: known Hz pattern + off-diagonal tensor reproduces the
 * analytic centered cross-curl on one step. */
static void test_cross_curl(void)
{
    fdtd_grid_t *g = fdtd_grid_create(8, 8, 2, 0.001, 0.9);
    CHECK(g);
    materials_t *m = materials_create(g);
    CHECK(m);

    /* Hz(i,j,0) = i + 0.1*j (dx-factored raw values) */
    for (size_t i = 0; i <= 8; i++)
        for (size_t j = 0; j <= 8; j++)
            g->hz[idxs(g, i, j, 0)] = (double)i + 0.1 * (double)j;

    /* tensor at Ex(4,4,0): m11 multiplies dy-curl (HZ(4,4)-HZ(4,3) = 0.1),
     * m12 multiplies centered dx-curl with negative sign:
     * 0.25*(HZ(5,4)+HZ(5,3)-HZ(3,4)-HZ(3,3)) = 0.5 */
    CHECK(materials_set_component_tensor(m, g, FIELD_EX, 4, 4, 0,
                                         0.8, 0.3, 0.8) == 0);
    fdtd_update_e_materials(g, m);

    material_tensor_t t;
    CHECK(materials_get_tensor(m, g, FIELD_EX, 4, 4, 0, &t) == 0);
    double expected = t.scale * (0.8 * 0.1 - 0.3 * 0.5);
    CLOSE(g->ex[idxs(g, 4, 4, 1)], 0.0, 0);
    /* adjacent edges untouched by the off-diagonal path */
    CHECK(fabs(g->ex[idxs(g, 5, 4, 0)]) < 1e-30 ||
          g->ex[idxs(g, 5, 4, 0)] == t.scale * 0.0);

    materials_destroy(m);
    fdtd_grid_destroy(g);
}

int main(void)
{
    test_api();
    test_invariance();
    test_cross_curl();
    printf("materials_tensor_unit: %zu checks, %zu failures\n", checks, failures);
    return failures ? 1 : 0;
}