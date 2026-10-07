#include <stdio.h>
#include <math.h>
#include "../engine/world/material_catalog.h"

#define CHECK(condition) do { if (!(condition)) { fprintf(stderr, "FAIL at line %d: %s\n", __LINE__, #condition); goto cleanup; } } while (0)

static int near(double a, double b)
{
    return isfinite(a) &&
           isfinite(b) &&
           fabs(a - b) <= 1e-12 * fmax(1.0, fmax(fabs(a), fabs(b)));
}

int main(void)
{
    int result = 1;
    material_catalog_t *cat = material_catalog_create_default();

    catalog_material_id_t id;
    catalog_material_id_t other;
    catalog_material_id_t old;
    catalog_material_id_t soil_id;

    catalog_properties_t p;
    catalog_properties_t q;
    catalog_definition_t d;

    CHECK(cat != NULL);
    CHECK(material_catalog_count(cat) == 23);

    /* Air and aliases. */
    CHECK(material_catalog_find(cat, "air", &id) == CATALOG_OK);
    CHECK(id == 0);

    CHECK(material_catalog_evaluate(
              cat, "vacuum", 915e6, 0, &p) == CATALOG_OK);
    CHECK(p.eps_r == 1.0);
    CHECK(p.sigma_s_per_m == 0.0);

    /* Generic concrete model and strict frequency range. */
    CHECK(material_catalog_evaluate(
              cat, "concrete_normal", 1e9, 0, &p) == CATALOG_OK);
    CHECK(near(p.eps_r, 5.24));
    CHECK(near(p.sigma_s_per_m, 0.0462));

    CHECK(material_catalog_evaluate(
              cat, "concrete", 2.4e9, 0, &q) == CATALOG_OK);
    CHECK(near(q.sigma_s_per_m,
               0.0462 * pow(2.4, 0.7822)));

    CHECK(material_catalog_evaluate(
              cat, "concrete", 915e6, 0, &q) ==
          CATALOG_FREQUENCY_RANGE);

    /* Proxies require permission. */
    CHECK(material_catalog_evaluate(
              cat, "wood_oak", 1e9, 0, &p) ==
          CATALOG_PROXY_DISABLED);

    CHECK(material_catalog_evaluate(
              cat, "wood_oak", 1e9, 1, &p) == CATALOG_OK);
    CHECK(p.status == CATALOG_PROXY);

    CHECK(material_catalog_evaluate(cat, "asphalt", 1e9, 1, &p) ==CATALOG_NEEDS_PROPERTIES);

    CHECK(material_catalog_evaluate(cat, "human_phantom", 1e9, 1, &p) == CATALOG_NEEDS_PROPERTIES);

    CHECK(material_catalog_find(cat, "dirt", &id) == CATALOG_OK);
    CHECK(material_catalog_find(cat, "soil_medium_dry", &other) == CATALOG_OK);
    CHECK(id == other);
    CHECK(material_catalog_find(cat, "soil_medium_dry", &soil_id) == CATALOG_OK);

    CHECK(material_catalog_evaluate(cat, "soil_wet", 1e9, 0, &p) == CATALOG_OK);
    CHECK(near(p.eps_r, 30.290) && near(p.sigma_s_per_m, 0.17151));
    CHECK(p.status == CATALOG_SOURCED);
    CHECK(material_catalog_evaluate(cat, "soil_medium_dry", 0.5e9, 0, &q) == CATALOG_FREQUENCY_RANGE);
    CHECK(material_catalog_evaluate(cat, "dirt", 2.4e9, 0, &q) == CATALOG_OK);
    CHECK(q.id == soil_id);

    CHECK(material_catalog_define_constant(cat, "my_wall", 4.0, 0.01, 1e9, "Test-only definition", "Unit test; not measured", 0, &id) == CATALOG_OK);

    CHECK(material_catalog_evaluate(cat, "my_wall", 1e9, 0, &p) == CATALOG_OK);
    CHECK(p.eps_r == 4.0);
    CHECK(p.status == CATALOG_USER_DEFINED);

    CHECK(material_catalog_define_constant(cat, "my_wall", 5.0, 0.02, 1e9, "Test", "Unit test", 0, &other) == CATALOG_EXISTS);

    CHECK(material_catalog_define_constant(cat, "my_wall", 5.0, 0.02, 1e9, "Test", "Unit test", 1, &other) == CATALOG_OK);
    CHECK(id == other);

    CHECK(material_catalog_find(cat, "asphalt", &old) == CATALOG_OK);

    CHECK(material_catalog_define_constant(cat, "asphalt", 4.0, 0.01, 1e9, "Illustration; not validated asphalt", "Unit test", 1, &id) == CATALOG_OK);
    CHECK(id == old);

    CHECK(material_catalog_get(cat, id, &d) == CATALOG_OK);
    CHECK(d.status == CATALOG_USER_DEFINED);

    CHECK(material_catalog_define_constant(cat, "air", 2.0, 0.0, 1e9, "Test", "Test", 1, &id) == CATALOG_RESERVED);

    CHECK(material_catalog_define_constant(cat, "dirt", 4.0, 0.0, 1e9, "Test", "Test", 1, &id) == CATALOG_RESERVED);

    CHECK(material_catalog_define_constant(cat, "bad", NAN, 0.0, 1e9, "Test", "Test", 0, &id) == CATALOG_INVALID);

    CHECK(material_catalog_define_constant(cat, "bad", 2.0, -1.0, 1e9, "Test", "Test", 0, &id) == CATALOG_INVALID);

    CHECK(material_catalog_evaluate(cat, "air", INFINITY, 0, &p) == CATALOG_INVALID);

    CHECK(material_catalog_find(cat, "unknown", &id) == CATALOG_NOT_FOUND);

    CHECK(material_catalog_find(NULL, "air", &id) == CATALOG_INVALID);

    CHECK(material_catalog_evaluate(cat, "air", 1e9, 0, NULL) == CATALOG_INVALID);

    while (material_catalog_count(cat) < MATERIAL_CATALOG_LIMIT) {
        char name[64];

        snprintf(name, sizeof(name), "test_%zu", material_catalog_count(cat));

        CHECK(material_catalog_define_constant(cat, name, 2.0, 0.0, 1e9, "Test", "Test", 0, &id) == CATALOG_OK);
    }

    CHECK(id == 255);

    CHECK(material_catalog_define_constant(cat, "overflow", 2.0, 0.0, 1e9, "Test", "Test", 0, &id) == CATALOG_FULL);

    printf("material_catalog_unit: PASS\n");
    result = 0;

cleanup:
    material_catalog_destroy(cat);
    return result;
}