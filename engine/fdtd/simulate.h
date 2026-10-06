#ifndef ENGINE_FDTD_SIMULATE_H
#define ENGINE_FDTD_SIMULATE_H

#include <stddef.h>
#include "field_component.h"

typedef struct {
    size_t nx, ny, nz;
    double dx;
    double courant_safety;
    size_t n_steps;

    size_t npml;

    // source: gaussian derivative, position already resolved to real cells */
    size_t src_i, src_j, src_k;
    field_component_t src_component;
    double t0_over_dt, tau_over_dt;

    size_t probe_i, probe_j, probe_k;
    field_component_t probe_component;
} fdtd_sim_config_t;

int fdtd_simulate_echo(const fdtd_sim_config_t *cfg, double *out_td, double *sample_rate_out);

#endif