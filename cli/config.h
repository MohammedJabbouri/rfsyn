#ifndef ENGINE_CLI_CONFIG_H
#define ENGINE_CLI_CONFIG_H

#include <stdio.h>

typedef struct config config_t;

config_t *config_create_empty(void);
config_t *config_load(const char *path);
int config_save(const config_t *cfg, const char *path);
void config_print(const config_t *cfg, FILE *out);
void config_destroy(config_t *cfg);

const char *config_get(const config_t *cfg, const char *section, const char *key);
double config_get_double(const config_t *cfg, const char *section, const char *key, double default_value);
long config_get_long(const config_t *cfg, const char *section, const char *key, long default_value);
int config_get_bool(const config_t *cfg, const char *section, const char *key, int default_value);

void config_set_string(config_t *cfg, const char *section, const char *key, const char *value);
void config_set_double(config_t *cfg, const char *section, const char *key, double value);
void config_set_long(config_t *cfg, const char *section, const char *key, long value);
void config_set_bool(config_t *cfg, const char *section, const char *key, int value);

#endif