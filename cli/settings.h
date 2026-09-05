#ifndef ENGINE_CLI_SETTINGS_H
#define ENGINE_CLI_SETTINGS_H

#include <stddef.h>
#include <stdint.h>
#include "config.h"
#include "../engine/core/transform.h"

int settings_build_chain(const config_t *cfg, uint64_t example_index, transform_t *stages_out, size_t max_stages, size_t *n_stages_out);

#endif