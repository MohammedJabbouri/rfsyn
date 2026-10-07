#include "material_config.h"
#include "cJSON.h"
#include <math.h>
#include <string.h>
#include <stdint.h>

#define MC_OK 0
#define MC_INVALID -1
#define MC_PARSE -2

static void log_line(FILE *log, const char *message) {
    if (log && message) fprintf(log, "%s\n", message);
}

static int copy_catalog_text(char *dst, size_t capacity, const char *src) {
    if (!src) src = "";
    size_t n = strlen(src);
    if (n >= capacity) return CATALOG_INVALID;
    memcpy(dst, src, n + 1);
    return CATALOG_OK;
}

static int parse_color(const char *text, color_rgb_t *out) {
    if (!out || !text || strlen(text) != 7 || text[0] != '#') return 0;
    unsigned value = 0;
    for (int i = 1; i <= 6; i++) {
        unsigned d;
        char ch = text[i];
        if (ch >= '0' && ch <= '9') d = (unsigned)(ch - '0');
        else if (ch >= 'A' && ch <= 'F') d = (unsigned)(ch - 'A') + 10u;
        else if (ch >= 'a' && ch <= 'f') d = (unsigned)(ch - 'a') + 10u;
        else return 0;
        value = (value << 4) | d;
    }
    out->r = (uint8_t)(value >> 16);
    out->g = (uint8_t)(value >> 8);
    out->b = (uint8_t)value;
    return 1;
}

static int parse_material_one(const char *name, const cJSON *node, material_catalog_t *catalog, FILE *log) {
    if (!node || !cJSON_IsObject(node)) {
        if (log) fprintf(log, "material_config: materials.%s must be object\n", name ? name : "(null)");
        return MC_INVALID;
    }
    catalog_definition_t d = {0};
    const cJSON *model = cJSON_GetObjectItemCaseSensitive((cJSON*)node, "model");
    const cJSON *eps = cJSON_GetObjectItemCaseSensitive((cJSON*)node, "eps_r");
    const cJSON *sig = cJSON_GetObjectItemCaseSensitive((cJSON*)node, "sigma_s_per_m");
    const cJSON *f0 = cJSON_GetObjectItemCaseSensitive((cJSON*)node, "reference_frequency_hz");
    const cJSON *descr = cJSON_GetObjectItemCaseSensitive((cJSON*)node, "description");
    const cJSON *src = cJSON_GetObjectItemCaseSensitive((cJSON*)node, "source");
    const cJSON *override_flag = cJSON_GetObjectItemCaseSensitive((cJSON*)node, "override");

    if (!cJSON_IsString(model) || strcmp(model->valuestring, "constant") != 0 || !cJSON_IsNumber(eps) || !cJSON_IsNumber(sig)) {
        if (log) fprintf(log, "material_config: %s requires model='constant' and numeric eps_r/sigma_s_per_m\n", name);
        return MC_INVALID;
    }
    if (copy_catalog_text(d.name, sizeof(d.name), name) != CATALOG_OK) return MC_INVALID;
    if (copy_catalog_text(d.description, sizeof(d.description), cJSON_IsString(descr) ? descr->valuestring : "") != CATALOG_OK) return MC_INVALID;
    if (copy_catalog_text(d.source, sizeof(d.source), cJSON_IsString(src) ? src->valuestring : "User-supplied values") != CATALOG_OK) return MC_INVALID;
    d.model = CATALOG_CONSTANT;
    d.status = CATALOG_USER_DEFINED;
    d.eps_r = eps->valuedouble;
    d.sigma_s_per_m = sig->valuedouble;
    d.reference_frequency_hz = cJSON_IsNumber(f0) ? f0->valuedouble : 1e9;
    return material_catalog_define(catalog, &d, cJSON_IsBool(override_flag) && cJSON_IsTrue(override_flag) ? 1 : 0, NULL);
}

int material_config_load(const cJSON *root, material_catalog_t *catalog, material_config_result_t *result, FILE *log) {
    if (!catalog || !result) return MC_INVALID;
    memset(result, 0, sizeof(*result));
    result->evaluation_frequency_hz = 1e9;
    if (!root) return MC_OK;
    if (!cJSON_IsObject(root)) {
        log_line(log, "material_config: root must be a JSON object");
        return MC_PARSE;
    }
    const cJSON *catalog_node = cJSON_GetObjectItemCaseSensitive((cJSON*)root, "material_catalog");
    if (catalog_node) {
        if (!cJSON_IsObject(catalog_node)) {
            log_line(log, "material_config: material_catalog must be object");
            return MC_INVALID;
        }
        const cJSON *enabled = cJSON_GetObjectItemCaseSensitive(catalog_node, "enabled");
        const cJSON *frequency = cJSON_GetObjectItemCaseSensitive(catalog_node, "evaluation_frequency_hz");
        const cJSON *allow = cJSON_GetObjectItemCaseSensitive(catalog_node, "allow_proxy");
        result->enabled = cJSON_IsBool(enabled) ? cJSON_IsTrue(enabled) : 0;
        result->allow_proxy = cJSON_IsBool(allow) ? cJSON_IsTrue(allow) : 0;
        if (cJSON_IsNumber(frequency)) result->evaluation_frequency_hz = frequency->valuedouble;
        if (!(result->evaluation_frequency_hz > 0.0 && isfinite(result->evaluation_frequency_hz))) {
            log_line(log, "material_config: evaluation_frequency_hz must be positive and finite");
            return MC_INVALID;
        }
    }
    const cJSON *materials_node = cJSON_GetObjectItemCaseSensitive((cJSON*)root, "materials");
    if (materials_node) {
        if (!cJSON_IsObject(materials_node)) {
            log_line(log, "material_config: materials must be object");
            return MC_INVALID;
        }
        for (const cJSON *entry = materials_node->child; entry; entry = entry->next) {
            if (!entry->string || !*entry->string) {
                log_line(log, "material_config: materials entry has empty key");
                return MC_INVALID;
            }
            int rc = parse_material_one(entry->string, entry, catalog, log);
            if (rc != CATALOG_OK) {
                if (log) fprintf(log, "material_config: %s failed: %s\n", entry->string, material_catalog_error(rc));
                return MC_INVALID;
            }
            result->materials_loaded++;
        }
    }
    if (result->enabled) {
        const cJSON *colors_node = cJSON_GetObjectItemCaseSensitive((cJSON*)root, "colors");
        if (colors_node) {
            if (!cJSON_IsObject(colors_node)) {
                log_line(log, "material_config: colors must be object");
                return MC_INVALID;
            }
            for (const cJSON *entry = colors_node->child; entry; entry = entry->next) {
                color_rgb_t rgb;
                if (!entry->string || !cJSON_IsString(entry)) {
                    log_line(log, "material_config: each color requires '#RRGGBB': 'material'");
                    return MC_INVALID;
                }
                if (!parse_color(entry->string, &rgb)) {
                    if (log) fprintf(log, "material_config: invalid color key '%s'\n", entry->string);
                    return MC_INVALID;
                }
                int rc = material_catalog_map_color(catalog, rgb, entry->valuestring);
                if (rc != CATALOG_OK) {
                    if (log) fprintf(log, "material_config: color %s -> %s failed: %s\n", entry->string, entry->valuestring, material_catalog_error(rc));
                    return MC_INVALID;
                }
                result->colors_seen++;
            }
        }
    }
    return MC_OK;
}

const char *material_config_error(int result) {
    switch (result) {
        case MC_OK: return "success";
        case MC_INVALID: return "invalid material configuration";
        case MC_PARSE: return "invalid JSON root";
        default: return "unknown material_config error";
    }
}
