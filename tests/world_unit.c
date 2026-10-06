#include <stdio.h>
#include <math.h>
#include "../engine/world/world.h"

#define CHECK(condition) do { if (!(condition)) { fprintf(stderr, "FAIL: line %d: %s\n", __LINE__, #condition); goto cleanup; } } while (0)

static int near(double a, double b) {
    return isfinite(a) && isfinite(b) && fabs(a - b) < 1e-12;
}

int main(void) {
    int result = 1;
    world_t *world = NULL;
    double heights[6] = {1.25, 1.25, 1.25, 1.25, 1.25, 1.25};
    world_material_id_t ids[6] = {1, 1, 1, 1, 1, 1};
    world_heightmap_config_t cfg = {3, 2, 4, 1.0, heights, ids, 0};
    world_cell_t cell;
    world = world_create_heightmap(&cfg);
    CHECK(world != NULL);
    CHECK(near(world_solid_volume_m3(world), 7.5));
    CHECK(world_mixed_cell_count(world) == 6);
    CHECK(world_get_cell(world, 0, 0, 0, &cell) == 0 && cell.count == 1 && cell.ids[0] == 1);
    CHECK(world_get_cell(world, 0, 0, 1, &cell) == 0 && cell.count == 2);
    CHECK(cell.ids[0] == 1 && cell.ids[1] == 0);
    CHECK(near(cell.fractions[0], 0.25) && near(cell.fractions[1], 0.75));
    CHECK(near(cell.normal[0], 0.0) && near(cell.normal[1], 0.0) && near(cell.normal[2], 1.0));
    CHECK(world_get_cell(world, 0, 0, 2, &cell) == 0 && cell.count == 1 && cell.ids[0] == 0);
    CHECK(world_get_cell(world, 3, 0, 0, &cell) == -1);
    CHECK(world_get_cell(world, 0, 0, 0, NULL) == -1);
    world_destroy(world); world = NULL;

    /* Samples of a linear ramp; reconstruction remains columnwise constant. */
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 2; j++) heights[i * 2 + j] = 0.25 + (double)i;
    }
    ids[2] = 2;
    world = world_create_heightmap(&cfg);
    CHECK(world != NULL && near(world_solid_volume_m3(world), 7.5));
    CHECK(world_get_cell(world, 1, 0, 1, &cell) == 0 && cell.count == 2 && cell.ids[0] == 2);
    CHECK(near(cell.normal[0], -1.0 / sqrt(2.0)) && near(cell.normal[2], 1.0 / sqrt(2.0)));
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 2; j++) {
            for (size_t k = 0; k < 4; k++) {
                CHECK(world_get_cell(world, i, j, k, &cell) == 0);
                double sum = 0.0;
                for (size_t n = 0; n < cell.count; n++) {
                    CHECK(cell.fractions[n] >= 0.0 && cell.fractions[n] <= 1.0);
                    sum += cell.fractions[n];
                }
                CHECK(near(sum, 1.0));
            }
        }
    }
    world_destroy(world); world = NULL;

    for (size_t n = 0; n < 6; n++) heights[n] = 0.0;
    heights[0] = 4.0;
    world = world_create_heightmap(&cfg);
    CHECK(world != NULL && world_mixed_cell_count(world) == 0 && near(world_solid_volume_m3(world), 4.0));
    world_destroy(world); world = NULL;
    cfg.dx = 0.5;
    for (size_t n = 0; n < 6; n++) heights[n] = 0.625;
    world = world_create_heightmap(&cfg);
    CHECK(world != NULL && near(world_solid_volume_m3(world), 0.9375));
    world_destroy(world); world = NULL;
    heights[0] = NAN;
    CHECK(world_create_heightmap(&cfg) == NULL);
    heights[0] = -1.0;
    CHECK(world_create_heightmap(&cfg) == NULL);
    heights[0] = 2.1;
    CHECK(world_create_heightmap(&cfg) == NULL);
    heights[0] = 0.5;
    ids[0] = 0;
    CHECK(world_create_heightmap(&cfg) == NULL);
    CHECK(world_create_heightmap(NULL) == NULL);
    printf("world_unit: PASS\n");
    result = 0;
cleanup:
    world_destroy(world);
    return result;
}