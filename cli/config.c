#include "config.h"
#include "cJSON.h"
#include <stdlib.h>
#include <string.h>

struct config {
    cJSON *root;
};

config_t *config_create_empty(void) {
    config_t *cfg = malloc(sizeof(config_t));
    if (!cfg) return NULL;
    cfg->root = cJSON_CreateObject();
    if (!cfg->root) { free(cfg); return NULL; }
    return cfg;
}

config_t *config_load(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) return NULL;

    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (size < 0) { fclose(f); return NULL; }

    char *buf = malloc((size_t)size + 1);
    if (!buf) { fclose(f); return NULL; }
    size_t read_len = fread(buf, 1, (size_t)size, f);
    fclose(f);
    buf[read_len] = '\0';

    cJSON *root = cJSON_Parse(buf);
    free(buf);
    if (!root) {
        const char *err = cJSON_GetErrorPtr();
        fprintf(stderr, "config: JSON parse error near: %s\n", err ? err : "(unknown)");
        return NULL;
    }

    config_t *cfg = malloc(sizeof(config_t));
    if (!cfg) { cJSON_Delete(root); return NULL; }
    cfg->root = root;
    return cfg;
}

int config_save(const config_t *cfg, const char *path) {
    if (!cfg) return -1;
    char *text = cJSON_Print(cfg->root);
    if (!text) return -1;

    FILE *f = fopen(path, "w");
    if (!f) { free(text); return -1; }
    size_t len = strlen(text);
    size_t written = fwrite(text, 1, len, f);
    fclose(f);
    free(text);
    return (written == len) ? 0 : -1;
}

void config_print(const config_t *cfg, FILE *out) {
    if (!cfg) return;
    char *text = cJSON_Print(cfg->root);
    if (!text) return;
    fprintf(out, "%s\n", text);
    free(text);
}

void config_destroy(config_t *cfg) {
    if (!cfg) return;
    cJSON_Delete(cfg->root);
    free(cfg);
}

const char *config_get(const config_t *cfg, const char *section, const char *key) {
    if (!cfg) return NULL;
    cJSON *s = cJSON_GetObjectItemCaseSensitive(cfg->root, section);
    if (!s) return NULL;
    cJSON *item = cJSON_GetObjectItemCaseSensitive(s, key);
    if (!cJSON_IsString(item)) return NULL;
    return item->valuestring;
}

double config_get_double(const config_t *cfg, const char *section, const char *key, double default_value) {
    if (!cfg) return default_value;
    cJSON *s = cJSON_GetObjectItemCaseSensitive(cfg->root, section);
    if (!s) return default_value;
    cJSON *item = cJSON_GetObjectItemCaseSensitive(s, key);
    if (!cJSON_IsNumber(item)) return default_value;
    return item->valuedouble;
}

long config_get_long(const config_t *cfg, const char *section, const char *key, long default_value) {
    return (long)config_get_double(cfg, section, key, (double)default_value);
}

int config_get_bool(const config_t *cfg, const char *section, const char *key, int default_value) {
    if (!cfg) return default_value;
    cJSON *s = cJSON_GetObjectItemCaseSensitive(cfg->root, section);
    if (!s) return default_value;
    cJSON *item = cJSON_GetObjectItemCaseSensitive(s, key);
    if (!cJSON_IsBool(item)) return default_value;
    return cJSON_IsTrue(item) ? 1 : 0;
}

static cJSON *get_or_create_section(cJSON *root, const char *section) {
    cJSON *s = cJSON_GetObjectItemCaseSensitive(root, section);
    if (s) return s;
    return cJSON_AddObjectToObject(root, section);
}

static void set_item(cJSON *section, const char *key, cJSON *newitem) {
    if (cJSON_HasObjectItem(section, key)) {
        cJSON_ReplaceItemInObjectCaseSensitive(section, key, newitem);
    } else {
        cJSON_AddItemToObject(section, key, newitem);
    }
}

void config_set_string(config_t *cfg, const char *section, const char *key, const char *value) {
    if (!cfg) return;
    cJSON *s = get_or_create_section(cfg->root, section);
    set_item(s, key, cJSON_CreateString(value));
}

void config_set_double(config_t *cfg, const char *section, const char *key, double value) {
    if (!cfg) return;
    cJSON *s = get_or_create_section(cfg->root, section);
    set_item(s, key, cJSON_CreateNumber(value));
}

void config_set_long(config_t *cfg, const char *section, const char *key, long value) {
    config_set_double(cfg, section, key, (double)value);
}

void config_set_bool(config_t *cfg, const char *section, const char *key, int value) {
    if (!cfg) return;
    cJSON *s = get_or_create_section(cfg->root, section);
    set_item(s, key, cJSON_CreateBool(value));
}