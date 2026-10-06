#ifndef ENGINE_CLI_FDTD_STAGE_H
#define ENGINE_CLI_FDTD_STAGE_H

#include "../engine/fdtd/simulate.h"
#include "../engine/core/transform.h"

transform_t fdtd_stage_create(const fdtd_sim_config_t *cfg);

#endif