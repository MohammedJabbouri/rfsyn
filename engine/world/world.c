#include "world.h"
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

typedef struct {
    size_t idx;
    world_material_id_t solid_id;
    double solid_fraction;
    double normal[3];
} mixed_cell_t;

struct world {
    size_t nx, ny, nz, total;
    double dx, solid_volume;
    world_material_id_t air_id;
    world_material_id_t *bulk_ids;
    mixed_cell_t *mixed;
    size_t mixed_count;
};

static size_t index_of(const world_t *world, size_t i, size_t j, size_t k) {
    return (i * world->ny + j) * world->nz + k;
}

static double fraction_at(double height_cells, size_t k) {
    double fraction = height_cells - (double)k;
    if (fraction <= 0.0) return 0.0;
    if (fraction >= 1.0) return 1.0;
    return fraction;
}

static int normal_at(const world_heightmap_config_t *cfg, size_t i, size_t j, double normal[3]) {
    size_t il = i ? i - 1 : i;
    size_t ir = i + 1 < cfg->nx ? i + 1 : i;
    size_t jl = j ? j - 1 : j;
    size_t jr = j + 1 < cfg->ny ? j + 1 : j;
    double gx = 0.0, gy = 0.0;
    if (ir != il) gx = (cfg->heights_m[ir * cfg->ny + j] - cfg->heights_m[il * cfg->ny + j]) / ((double)(ir - il) * cfg->dx);
    if (jr != jl) gy = (cfg->heights_m[i * cfg->ny + jr] - cfg->heights_m[i * cfg->ny + jl]) / ((double)(jr - jl) * cfg->dx);
    double length = hypot(hypot(gx, gy), 1.0);
    if (!isfinite(length) || length <= 0.0) return -1;
    normal[0] = -gx / length;
    normal[1] = -gy / length;
    normal[2] = 1.0 / length;
    return 0;
}

world_t *world_create_heightmap(const world_heightmap_config_t *cfg) {
    if (!cfg || !cfg->heights_m || !cfg->column_materials) return NULL;
    if (!cfg->nx || !cfg->ny || !cfg->nz || !isfinite(cfg->dx) || cfg->dx <= 0.0) return NULL;
    if (cfg->nx > SIZE_MAX / cfg->ny) return NULL;
    size_t columns = cfg->nx * cfg->ny;
    if (columns > SIZE_MAX / cfg->nz) return NULL;
    size_t total = columns * cfg->nz;
    double domain_height = (double)cfg->nz * cfg->dx;
    double cell_volume = cfg->dx * cfg->dx * cfg->dx;
    if (!isfinite(domain_height) || !isfinite(cell_volume) || cell_volume <= 0.0) return NULL;

    size_t mixed_count = 0;
    for (size_t col = 0; col < columns; col++) {
        double h = cfg->heights_m[col];
        if (!isfinite(h) || h < 0.0 || h > domain_height) return NULL;
        if (h > 0.0 && cfg->column_materials[col] == cfg->air_id) return NULL;
        double hc = fmin(h / cfg->dx, (double)cfg->nz);
        if (hc > 0.0 && hc < (double)cfg->nz && hc != floor(hc)) mixed_count++;
    }
    if (mixed_count > SIZE_MAX / sizeof(mixed_cell_t)) return NULL;

    world_t *world = calloc(1, sizeof(*world));
    if (!world) return NULL;
    world->nx = cfg->nx; world->ny = cfg->ny; world->nz = cfg->nz;
    world->dx = cfg->dx; world->air_id = cfg->air_id; world->total = total;
    world->bulk_ids = malloc(total * sizeof(*world->bulk_ids));
    if (!world->bulk_ids) { world_destroy(world); return NULL; }
    if (mixed_count) {
        world->mixed = malloc(mixed_count * sizeof(*world->mixed));
        if (!world->mixed) { world_destroy(world); return NULL; }
    }

    for (size_t i = 0; i < cfg->nx; i++) {
        for (size_t j = 0; j < cfg->ny; j++) {
            size_t col = i * cfg->ny + j;
            double hc = fmin(cfg->heights_m[col] / cfg->dx, (double)cfg->nz);
            for (size_t k = 0; k < cfg->nz; k++) {
                size_t idx = index_of(world, i, j, k);
                double fraction = fraction_at(hc, k);
                world->bulk_ids[idx] = fraction >= 0.5 ? cfg->column_materials[col] : cfg->air_id;
                world->solid_volume += fraction * cell_volume;
                if (fraction > 0.0 && fraction < 1.0) {
                    mixed_cell_t *m = &world->mixed[world->mixed_count++];
                    m->idx = idx;
                    m->solid_id = cfg->column_materials[col];
                    m->solid_fraction = fraction;
                    if (normal_at(cfg, i, j, m->normal) != 0) { world_destroy(world); return NULL; }
                }
            }
        }
    }
    if (!isfinite(world->solid_volume)) { world_destroy(world); return NULL; }
    return world;
}

void world_destroy(world_t *world) {
    if (!world) return;
    free(world->bulk_ids);
    free(world->mixed);
    free(world);
}

int world_get_cell(const world_t *world, size_t i, size_t j, size_t k, world_cell_t *out) {
    if (!world || !out || i >= world->nx || j >= world->ny || k >= world->nz) return -1;
    size_t idx = index_of(world, i, j, k);
    *out = (world_cell_t){0};
    size_t lo = 0, hi = world->mixed_count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (world->mixed[mid].idx < idx) lo = mid + 1;
        else hi = mid;
    }
    if (lo < world->mixed_count && world->mixed[lo].idx == idx) {
        const mixed_cell_t *m = &world->mixed[lo];
        out->count = 2;
        out->ids[0] = m->solid_id;
        out->ids[1] = world->air_id;
        out->fractions[0] = m->solid_fraction;
        out->fractions[1] = 1.0 - m->solid_fraction;
        for (size_t n = 0; n < 3; n++) out->normal[n] = m->normal[n];
    } else {
        out->count = 1;
        out->ids[0] = world->bulk_ids[idx];
        out->fractions[0] = 1.0;
    }
    return 0;
}

size_t world_mixed_cell_count(const world_t *world) {
    return world ? world->mixed_count : 0;
}

double world_solid_volume_m3(const world_t *world) {
    return world ? world->solid_volume : 0.0;
}