#include <stdio.h>
#include <math.h>
#include "../engine/fdtd/materials.h"
#include "../engine/core/constants.h"

#define CHECK(condition) do { if (!(condition)) { fprintf(stderr, "FAIL: line %d: %s\n", __LINE__, #condition); goto cleanup; } } while (0)

static int near(double a, double b) {
    return isfinite(a) && isfinite(b) && fabs(a - b) <= 1e-12 * fmax(1.0, fmax(fabs(a), fabs(b)));
}

int main(void) {
    int result = 1;
    fdtd_grid_t *g = fdtd_grid_create(12, 12, 12, 1e-3, 0.99);
    materials_t *mat = NULL;
    cpml_t *pml = NULL;
    material_coeff_t vacuum, ex, ey, ez;
    material_id_t id, again;
    CHECK(g != NULL);
    mat = materials_create(g);
    CHECK(mat != NULL);
    CHECK(materials_get_coefficients(mat, g, FIELD_EX, 6, 6, 6, &vacuum) == 0);
    CHECK(vacuum.Ca == 1.0);
    CHECK(near(vacuum.Cb, g->dt / (EPS0 * g->dx)));
    CHECK(materials_register(mat, 4.0, 0.0, &id) == 0);
    CHECK(materials_register(mat, 4.0, 0.0, &again) == 0 && id == again);
    CHECK(materials_assign_box(mat, g, 5, 8, 5, 8, 5, 8, id) == 0);
    CHECK(materials_get_coefficients(mat, g, FIELD_EY, 6, 6, 6, &ey) == 0);
    CHECK(near(ey.Cb, vacuum.Cb / 4.0));
    CHECK(materials_get_coefficients(mat, g, FIELD_EX, 4, 6, 6, &ex) == 0);
    CHECK(near(ex.Cb, vacuum.Cb));
    CHECK(materials_set_component(mat, g, FIELD_EX, 6, 6, 6, 9.0, 0.1) == 0);
    CHECK(materials_set_component(mat, g, FIELD_EZ, 6, 6, 6, 16.0, 0.0) == 0);
    CHECK(materials_get_coefficients(mat, g, FIELD_EX, 6, 6, 6, &ex) == 0);
    CHECK(materials_get_coefficients(mat, g, FIELD_EY, 6, 6, 6, &ey) == 0);
    CHECK(materials_get_coefficients(mat, g, FIELD_EZ, 6, 6, 6, &ez) == 0);
    CHECK(ex.Ca < 1.0 && ey.Ca == 1.0 && ez.Ca == 1.0);
    CHECK(near(ey.Cb, vacuum.Cb / 4.0) && near(ez.Cb, vacuum.Cb / 16.0));
    double loss = 0.1 * g->dt / (2.0 * 9.0 * EPS0);
    CHECK(near(ex.Ca, (1.0 - loss) / (1.0 + loss)));
    CHECK(near(ex.Cb, (vacuum.Cb / 9.0) / (1.0 + loss)));
    size_t idx = fdtd_index(g, 6, 6, 6);
    g->ex[idx] = g->ey[idx] = g->ez[idx] = 1.0;
    fdtd_update_e_materials(g, mat);
    CHECK(near(g->ex[idx], ex.Ca) && g->ey[idx] == 1.0 && g->ez[idx] == 1.0);
    pml = cpml_create(g, cpml_default_params(2));
    CHECK(pml != NULL);
    g->ex[idx] = g->ey[idx] = g->ez[idx] = 1.0;
    fdtd_update_e_cpml_materials(g, pml, mat);
    CHECK(near(g->ex[idx], ex.Ca) && g->ey[idx] == 1.0 && g->ez[idx] == 1.0);
    CHECK(materials_clear_component(mat, g, FIELD_EX, 6, 6, 6) == 0);
    CHECK(materials_get_coefficients(mat, g, FIELD_EX, 6, 6, 6, &ex) == 0 && near(ex.Cb, ey.Cb));
    CHECK(materials_set_box(mat, g, 6, 7, 6, 7, 6, 7, 1.0, 0.0) == 0);
    CHECK(materials_get_coefficients(mat, g, FIELD_EZ, 6, 6, 6, &ez) == 0 && near(ez.Cb, vacuum.Cb));
    CHECK(materials_set_box(mat, g, 0, 14, 0, 1, 0, 1, 4.0, 0.0) == -1);
    CHECK(materials_register(mat, NAN, 0.0, &id) == -1);
    CHECK(materials_register(mat, 4.0, -1.0, &id) == -1);
    CHECK(materials_set_component(mat, g, FIELD_HX, 6, 6, 6, 4.0, 0.0) == -1);
    CHECK(materials_set_component(mat, g, FIELD_EX, 12, 6, 6, 4.0, 0.0) == -1);
    CHECK(materials_get_coefficients(mat, g, FIELD_EX, 6, 6, 6, NULL) == -1);
    printf("materials_unit: PASS\n");
    result = 0;
cleanup:
    cpml_destroy(pml);
    materials_destroy(mat);
    fdtd_grid_destroy(g);
    return result;
}
