#include "material_catalog.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MATERIAL_CATALOG_COLOR_LIMIT 256

typedef struct {
    color_rgb_t color;
    char name[MATERIAL_CATALOG_NAME_SIZE];
} color_entry_t;

struct material_catalog {
    catalog_definition_t entries[MATERIAL_CATALOG_LIMIT];
    size_t count;

    color_entry_t colors[MATERIAL_CATALOG_COLOR_LIMIT];
    size_t color_count;
};

typedef struct { const char *name, *target; } alias_t;
static const alias_t aliases[] = {
    {"vacuum", "air"}, {"drywall", "plasterboard"}, {"dirt", "soil_medium_dry"},
    {"grass_dry_proxy", "soil_very_dry"}, {"grass_moist_proxy", "soil_medium_dry"},
    {"human_phantom", "tissue_muscle"}, {"wood", "wood_framing_spf"},
    {"brick", "brick_clay"}, {"concrete", "concrete_normal"}
};

static const char *canonical(const char *name) {
    if (!name) return NULL;
    for (size_t i = 0; i < sizeof(aliases) / sizeof(aliases[0]); i++) {
        if (strcmp(name, aliases[i].name) == 0) return aliases[i].target;
    }
    return name;
}

static int copy_text(char *dst, size_t capacity, const char *src) {
    if (!src) src = "";
    size_t n = strlen(src);
    if (n >= capacity) return CATALOG_INVALID;
    memcpy(dst, src, n + 1);
    return CATALOG_OK;
}

static int name_valid(const char *name) {
    if (!name || !*name) return 0;
    for (const unsigned char *p = (const unsigned char *)name; *p; p++) {
        if (!((*p >= 'a' && *p <= 'z') || (*p >= '0' && *p <= '9') || *p == '_')) return 0;
    }
    return 1;
}

static int definition_valid(const catalog_definition_t *d) {
    if (!d || !memchr(d->name, 0, sizeof(d->name)) || !memchr(d->description, 0, sizeof(d->description)) || !memchr(d->source, 0, sizeof(d->source))) return 0;
    if (!name_valid(d->name) || !d->source[0]) return 0;
    if (d->status < CATALOG_ANALYTIC || d->status > CATALOG_PENDING) return 0;
    if (!isfinite(d->min_frequency_hz) || !isfinite(d->max_frequency_hz)) return 0;
    if (!(d->min_frequency_hz == 0.0 && d->max_frequency_hz == 0.0) && !(d->min_frequency_hz > 0.0 && d->max_frequency_hz >= d->min_frequency_hz)) return 0;
    if (d->model == CATALOG_UNCONFIGURED) return d->status == CATALOG_PENDING;
    if (d->status == CATALOG_PENDING) return 0;
    if (d->model == CATALOG_CONSTANT) return isfinite(d->eps_r) && d->eps_r >= 1.0 && isfinite(d->sigma_s_per_m) && d->sigma_s_per_m >= 0.0 && isfinite(d->reference_frequency_hz) && d->reference_frequency_hz > 0.0;
    if (d->model == CATALOG_POWER_LAW) return isfinite(d->a) && d->a > 0.0 && isfinite(d->b) && isfinite(d->c) && d->c >= 0.0 && isfinite(d->d) && d->min_frequency_hz > 0.0;
    return 0;
}

int material_catalog_find(const material_catalog_t *catalog, const char *name, catalog_material_id_t *id_out) {
    if (!catalog || !name || !id_out) return CATALOG_INVALID;
    name = canonical(name);
    for (size_t i = 0; i < catalog->count; i++) {
        if (strcmp(catalog->entries[i].name, name) == 0) { *id_out = (catalog_material_id_t)i; return CATALOG_OK; }
    }
    return CATALOG_NOT_FOUND;
}

int material_catalog_get(const material_catalog_t *catalog, catalog_material_id_t id, catalog_definition_t *out) {
    if (!catalog || !out) return CATALOG_INVALID;
    if ((size_t)id >= catalog->count) return CATALOG_NOT_FOUND;
    *out = catalog->entries[id];
    return CATALOG_OK;
}

int material_catalog_define(material_catalog_t *catalog, const catalog_definition_t *definition, int override_existing, catalog_material_id_t *id_out) {
    if (!catalog || !id_out || !definition_valid(definition)) return CATALOG_INVALID;
    if (strcmp(canonical(definition->name), definition->name) != 0) return CATALOG_RESERVED;
    catalog_material_id_t id;
    int found = material_catalog_find(catalog, definition->name, &id);
    if (found == CATALOG_OK) {
        if (strcmp(definition->name, "air") == 0) return CATALOG_RESERVED;
        if (!override_existing) return CATALOG_EXISTS;
        catalog->entries[id] = *definition;
    } else {
        if (catalog->count == MATERIAL_CATALOG_LIMIT) return CATALOG_FULL;
        id = (catalog_material_id_t)catalog->count;
        catalog->entries[catalog->count++] = *definition;
    }
    *id_out = id;
    return CATALOG_OK;
}

int material_catalog_define_constant(material_catalog_t *catalog, const char *name, double eps_r, double sigma_s_per_m, double reference_frequency_hz, const char *description, const char *source, int override_existing, catalog_material_id_t *id_out) {
    catalog_definition_t d = {0};
    if (!name || copy_text(d.name, sizeof(d.name), name) != CATALOG_OK || copy_text(d.description, sizeof(d.description), description) != CATALOG_OK || copy_text(d.source, sizeof(d.source), source ? source : "User supplied; unverified") != CATALOG_OK) return CATALOG_INVALID;
    d.model = CATALOG_CONSTANT;
    d.status = CATALOG_USER_DEFINED;
    d.eps_r = eps_r;
    d.sigma_s_per_m = sigma_s_per_m;
    d.reference_frequency_hz = reference_frequency_hz;
    return material_catalog_define(catalog, &d, override_existing, id_out);
}

int material_catalog_evaluate(const material_catalog_t *catalog, const char *name, double frequency_hz, int allow_proxy, catalog_properties_t *out) {
    if (!out || !isfinite(frequency_hz) || frequency_hz <= 0.0) return CATALOG_INVALID;
    catalog_material_id_t id;
    int rc = material_catalog_find(catalog, name, &id);
    if (rc != CATALOG_OK) return rc;
    const catalog_definition_t *d = &catalog->entries[id];
    if (d->model == CATALOG_UNCONFIGURED) return CATALOG_NEEDS_PROPERTIES;
    if (d->status == CATALOG_PROXY && !allow_proxy) return CATALOG_PROXY_DISABLED;
    if (d->min_frequency_hz > 0.0 && (frequency_hz < d->min_frequency_hz || frequency_hz > d->max_frequency_hz)) return CATALOG_FREQUENCY_RANGE;
    double eps = d->eps_r, sigma = d->sigma_s_per_m;
    if (d->model == CATALOG_POWER_LAW) {
        double f = frequency_hz / 1e9;
        eps = d->a * pow(f, d->b);
        sigma = d->c == 0.0 ? 0.0 : d->c * pow(f, d->d);
    }
    if (!isfinite(eps) || eps < 1.0 || !isfinite(sigma) || sigma < 0.0) return CATALOG_INVALID;
    *out = (catalog_properties_t){id, d->status, frequency_hz, eps, sigma};
    return CATALOG_OK;
}

static int add_power(material_catalog_t *cat, const char *name, const char *description, catalog_status_t status, double a, double c, double exponent, double low, double high) {
    catalog_definition_t d = {0};
    catalog_material_id_t id;
    if (copy_text(d.name, sizeof(d.name), name) || copy_text(d.description, sizeof(d.description), description)) return -1;
    if (copy_text(d.source, sizeof(d.source), "ITU-R P.2040-3 Table 3; generic class model (MathWorks reproduction)")) return -1;
    d.model = CATALOG_POWER_LAW; d.status = status;
    d.a = a; d.b = 0.0; d.c = c; d.d = exponent;
    d.min_frequency_hz = low; d.max_frequency_hz = high;
    return material_catalog_define(cat, &d, 0, &id);
}

static int add_const(material_catalog_t *cat, const char *name, const char *description, catalog_status_t status, const char *source, double eps_r, double sigma, double ref_hz, double low_hz, double high_hz) {
    catalog_definition_t d = {0};
    catalog_material_id_t id;
    if (copy_text(d.name, sizeof(d.name), name) ||
        copy_text(d.description, sizeof(d.description), description) ||
        copy_text(d.source, sizeof(d.source), source)) return -1;
    d.model = CATALOG_CONSTANT; d.status = status;
    d.eps_r = eps_r; d.sigma_s_per_m = sigma; d.reference_frequency_hz = ref_hz;
    d.min_frequency_hz = low_hz; d.max_frequency_hz = high_hz;
    return material_catalog_define(cat, &d, 0, &id);
}

static int add_pending(material_catalog_t *cat, const char *name, const char *description) {
    catalog_definition_t d = {0};
    catalog_material_id_t id;
    if (copy_text(d.name, sizeof(d.name), name) || copy_text(d.description, sizeof(d.description), description)) return -1;
    if (copy_text(d.source, sizeof(d.source), "Not configured: supply documented properties before simulation")) return -1;
    d.model = CATALOG_UNCONFIGURED; d.status = CATALOG_PENDING;
    return material_catalog_define(cat, &d, 0, &id);
}

material_catalog_t *material_catalog_create_default(void) {
    material_catalog_t *cat = calloc(1, sizeof(*cat));
    if (!cat) return NULL;
    catalog_definition_t air = {0};
    catalog_material_id_t id;
    copy_text(air.name, sizeof(air.name), "air");
    copy_text(air.description, sizeof(air.description), "Air represented as vacuum; nonmagnetic");
    copy_text(air.source, sizeof(air.source), "Analytic vacuum approximation");
    air.model = CATALOG_CONSTANT; air.status = CATALOG_ANALYTIC;
    air.eps_r = 1.0; air.reference_frequency_hz = 1e9;
    if (material_catalog_define(cat, &air, 0, &id) != CATALOG_OK) goto fail;
    if (add_power(cat, "wood_framing_spf", "Generic wood proxy; not SPF-specific, isotropic", CATALOG_PROXY, 1.99, 0.0047, 1.0718, 1e6, 100e9)) goto fail;
    if (add_power(cat, "wood_pine", "Generic wood proxy; not pine-specific or a living-tree model", CATALOG_PROXY, 1.99, 0.0047, 1.0718, 1e6, 100e9)) goto fail;
    if (add_power(cat, "wood_oak", "Generic wood proxy; not oak-specific", CATALOG_PROXY, 1.99, 0.0047, 1.0718, 1e6, 100e9)) goto fail;
    if (add_power(cat, "brick_clay", "Generic brick proxy; not a clay-product-specific measurement", CATALOG_PROXY, 3.91, 0.0238, 0.16, 1e9, 10e9)) goto fail;
    if (add_power(cat, "brick_concrete", "Generic concrete proxy; not a masonry-brick-specific measurement", CATALOG_PROXY, 5.24, 0.0462, 0.7822, 1e9, 100e9)) goto fail;
    if (add_power(cat, "concrete_normal", "Generic concrete class; moisture and aggregate not specified", CATALOG_SOURCED, 5.24, 0.0462, 0.7822, 1e9, 100e9)) goto fail;
    if (add_pending(cat, "concrete_lightweight", "Lightweight concrete; product-specific properties required")) goto fail;
    if (add_pending(cat, "asphalt", "Asphalt pavement; mix and moisture assumptions required")) goto fail;
    if (add_pending(cat, "ceramic_tile", "Ceramic tile; product-specific properties required")) goto fail;
    if (add_pending(cat, "vinyl_flooring", "Vinyl covering; composition-specific properties required")) goto fail;
    if (add_pending(cat, "carpet", "Carpet layer; fibers, backing and effective properties required")) goto fail;
    if (add_pending(cat, "gravel", "Aggregate/air/water effective material; composition required")) goto fail;
    if (add_const(cat, "soil_very_dry", "Silty loam, mv=0.07, T=23 C; 1 GHz point value (sigma varies ~2x over 0.9-2.45 GHz)", CATALOG_SOURCED,  "ITU-R P.527-6 Sec. 5.2, eqs (57)-(70), Table 2 texture",  4.280, 0.02665, 1e9, 0.8e9, 2.6e9)) goto fail;
    if (add_const(cat, "soil_medium_dry",  "Silty loam, mv=0.20, T=23 C; 1 GHz point value (sigma varies ~2x over 0.9-2.45 GHz)", CATALOG_SOURCED, "ITU-R P.527-6 Sec. 5.2, eqs (57)-(70), Table 2 texture", 9.896, 0.06674, 1e9, 0.8e9, 2.6e9)) goto fail;
    if (add_const(cat, "soil_wet", "Silty loam, mv=0.50, T=23 C; 1 GHz point value (sigma varies ~2x over 0.9-2.45 GHz)", CATALOG_SOURCED, "ITU-R P.527-6 Sec. 5.2, eqs (57)-(70), Table 2 texture", 30.290, 0.17151, 1e9, 0.8e9, 2.6e9)) goto fail;
    if (add_power(cat, "plasterboard", "Generic plasterboard class", CATALOG_SOURCED, 2.73, 0.0085, 0.9395, 1e9, 100e9)) goto fail;
    if (add_power(cat, "glass", "Generic glass class; lower-frequency branch only", CATALOG_SOURCED, 6.31, 0.0036, 1.3394, 0.1e9, 100e9)) goto fail;
    if (add_power(cat, "ceiling_board", "Generic ceiling-board class; lower-frequency branch only", CATALOG_SOURCED, 1.48, 0.0011, 1.0750, 1e9, 100e9)) goto fail;

    if (add_pending(cat, "tissue_muscle", "Frequency-specific published tissue properties required"))
        goto fail;

    if (add_pending(cat, "tissue_fat", "Frequency-specific published tissue properties required"))
        goto fail;

    if (add_pending(cat, "tissue_skin", "Frequency-specific published tissue properties required"))
        goto fail;

    return cat;

fail:
    material_catalog_destroy(cat);
    return NULL;
}

void material_catalog_destroy(material_catalog_t *catalog)
{
    free(catalog);
}

size_t material_catalog_count(const material_catalog_t *catalog)
{
    return catalog ? catalog->count : 0;
}

int material_catalog_map_color(material_catalog_t *catalog,
                               color_rgb_t color,
                               const char *material_name)
{
    size_t slot;
    catalog_material_id_t id;

    if (!catalog ||
        !material_name ||
        strlen(material_name) >= MATERIAL_CATALOG_NAME_SIZE) {
        return CATALOG_INVALID;
    }

    if (material_catalog_find(catalog, material_name, &id) != CATALOG_OK) {
        return CATALOG_NOT_FOUND;
    }

    for (slot = 0; slot < catalog->color_count; slot++) {
        color_rgb_t stored = catalog->colors[slot].color;

        if (stored.r == color.r &&
            stored.g == color.g &&
            stored.b == color.b) {
            strcpy(catalog->colors[slot].name,
                   catalog->entries[id].name);
            return CATALOG_OK;
        }
    }

    if (catalog->color_count == MATERIAL_CATALOG_COLOR_LIMIT) {
        return CATALOG_FULL;
    }

    slot = catalog->color_count++;
    catalog->colors[slot].color = color;
    strcpy(catalog->colors[slot].name, catalog->entries[id].name);

    return CATALOG_OK;
}

int material_catalog_lookup_color(const material_catalog_t *catalog, color_rgb_t color, char *name_out, size_t name_capacity) {
    if (!catalog || !name_out || name_capacity == 0) {
        return CATALOG_INVALID;
    }

    for (size_t slot = 0; slot < catalog->color_count; slot++) {
        color_rgb_t stored = catalog->colors[slot].color;

        if (stored.r == color.r &&
            stored.g == color.g &&
            stored.b == color.b) {
            if (strlen(catalog->colors[slot].name) >= name_capacity) {
                return CATALOG_INVALID;
            }

            strcpy(name_out, catalog->colors[slot].name);
            return CATALOG_OK;
        }
    }

    return CATALOG_NOT_FOUND;
}

const char *material_catalog_error(int result)
{
    switch (result) {
        case CATALOG_OK:
            return "success";

        case CATALOG_INVALID:
            return "invalid definition, frequency, or argument";

        case CATALOG_NOT_FOUND:
            return "unknown material";

        case CATALOG_EXISTS:
            return "material exists; explicit override required";

        case CATALOG_FULL:
            return "catalog limit of 256 entries reached";

        case CATALOG_NEEDS_PROPERTIES:
            return "material has no configured physical properties";

        case CATALOG_PROXY_DISABLED:
            return "generic proxy requires explicit permission";

        case CATALOG_FREQUENCY_RANGE:
            return "frequency outside configured model range";

        case CATALOG_RESERVED:
            return "reserved air or alias name; use canonical material name";

        default:
            return "unknown catalog error";
    }
}