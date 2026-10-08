#ifndef ENGINE_CLI_MATERIAL_CONFIG_H
#define ENGINE_CLI_MATERIAL_CONFIG_H

#include <stdio.h>
#include "../engine/world/material_catalog.h"

typedef struct cJSON cJSON;

typedef struct {
    int enabled;
    int allow_proxy;
    double evaluation_frequency_hz;
    int materials_loaded;
    int colors_seen;
} material_config_result_t;

int material_config_load(const cJSON *root, material_catalog_t *catalog, material_config_result_t *result, FILE *log);

// LOADING IS NOT TRANSACTIONAL, ON FAILURE DISCARD THE CATALOG AND DO NOT RESULT, CREATE FRESH CATALOG

const char *material_config_error(int result);

#endif