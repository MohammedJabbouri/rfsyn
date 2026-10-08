#ifndef ENGINE_WORLD_MATERIAL_CATALOG_H
#define ENGINE_WORLD_MATERIAL_CATALOG_H
#include <stddef.h>
#include <stdint.h>
#define MATERIAL_CATALOG_LIMIT 256
#define MATERIAL_CATALOG_COLOR_LIMIT 256
#define MATERIAL_CATALOG_NAME_SIZE 64
#define MATERIAL_CATALOG_TEXT_SIZE 192

typedef struct material_catalog material_catalog_t;
typedef uint8_t catalog_material_id_t;
typedef enum { CATALOG_UNCONFIGURED, CATALOG_CONSTANT, CATALOG_POWER_LAW } catalog_model_t;
typedef enum { CATALOG_ANALYTIC, CATALOG_SOURCED, CATALOG_PROXY, CATALOG_USER_DEFINED, CATALOG_PENDING } catalog_status_t;
typedef enum { CATALOG_OK=0, CATALOG_INVALID=-1, CATALOG_NOT_FOUND=-2, CATALOG_EXISTS=-3, CATALOG_FULL=-4, CATALOG_NEEDS_PROPERTIES=-5, CATALOG_PROXY_DISABLED=-6, CATALOG_FREQUENCY_RANGE=-7, CATALOG_RESERVED=-8 } catalog_result_t;
typedef struct {
    char name[MATERIAL_CATALOG_NAME_SIZE];
    char description[MATERIAL_CATALOG_TEXT_SIZE];
    char source[MATERIAL_CATALOG_TEXT_SIZE];
    catalog_model_t model;
    catalog_status_t status;
    double eps_r, sigma_s_per_m, reference_frequency_hz;
    double a, b, c, d;
    double min_frequency_hz, max_frequency_hz;
} catalog_definition_t;
typedef struct {
    catalog_material_id_t id;
    catalog_status_t status;
    double frequency_hz, eps_r, sigma_s_per_m;
} catalog_properties_t;
typedef struct { uint8_t r, g, b; } color_rgb_t;
material_catalog_t *material_catalog_create_default(void);
void material_catalog_destroy(material_catalog_t *catalog);
size_t material_catalog_count(const material_catalog_t *catalog);
int material_catalog_find(const material_catalog_t *catalog, const char *name, catalog_material_id_t *id_out);
int material_catalog_get(const material_catalog_t *catalog, catalog_material_id_t id, catalog_definition_t *out);
int material_catalog_define(material_catalog_t *catalog, const catalog_definition_t *definition, int override_existing, catalog_material_id_t *id_out);
int material_catalog_define_constant(material_catalog_t *catalog, const char *name, double eps_r, double sigma_s_per_m, double reference_frequency_hz, const char *description, const char *source, int override_existing, catalog_material_id_t *id_out);
int material_catalog_evaluate(const material_catalog_t *catalog, const char *name, double frequency_hz, int allow_proxy, catalog_properties_t *out);
int material_catalog_map_color(material_catalog_t *catalog, color_rgb_t color, const char *material_name);
int material_catalog_lookup_color(const material_catalog_t *catalog, color_rgb_t color, char *name_out, size_t name_capacity);
const char *material_catalog_error(int result);
#endif