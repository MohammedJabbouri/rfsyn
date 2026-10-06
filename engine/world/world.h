#ifndef ENGINE_WORLD_WORLD_H
#define ENGINE_WORLD_WORLD_H

#include <stddef.h>
#include <stdint.h>

typedef uint8_t world_material_id_t;
typedef struct world world_t;

typedef struct {
    size_t nx, ny, nz;
    double dx;

    const double *heights_m;
    const world_material_id_t *column_materials;
    world_material_id_t air_id;
} world_heightmap_config_t;

typedef struct {
    world_material_id_t ids[2];
    double fractions[2];
    size_t count;
    double normal[3];
} world_cell_t;

world_t *world_create_heightmap(const world_heightmap_config_t *cfg);
void world_destroy(world_t *world);
int world_get_cell(const world_t *world, size_t i, size_t j, size_t k, world_cell_t *out);
size_t world_mixed_cell_count(const world_t *world);
double world_solid_volume_m3(const world_t *world);

#endif
