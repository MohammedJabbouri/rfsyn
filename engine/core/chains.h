#ifndef ENGINE_CORE_CHAINS_H
#define ENGINE_CORE_CHAINS_H

#include <stddef.h>
#include "rf_signal.h"
#include "transform.h"


typedef struct chains_registry chains_registry_t;

chains_registry_t *chains_registry_create(void);
void chains_registry_destroy(chains_registry_t *reg);

int chains_register(chains_registry_t *reg, const char *name, const transform_t *stages, size_t n_stages);

int chains_lookup(const chains_registry_t *reg, const char *name, const transform_t **stages_out, size_t *n_stages_out);

int chains_run(const chains_registry_t *reg, const char *name, signal_t *sig);

#endif