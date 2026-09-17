#ifndef ENGINE_FDTD_PROBE_H
#define ENGINE_FDTD_PROBE_H

#include <stddef.h>
#include "grid.h"
#include "field_component.h"

typedef struct {
    size_t i, j, k;
    field_component_t component;
    double *buffer;
    size_t n_steps;
    size_t next_index;
} probe_td_t;

int  probe_td_create(probe_td_t *probe, size_t i, size_t j, size_t k,
                      field_component_t component, size_t n_steps);
void probe_td_destroy(probe_td_t *probe);

void probe_td_record(probe_td_t *probe, const fdtd_grid_t *g);

typedef struct {
    size_t i, j, k;
    field_component_t component;

    const double *freqs_hz;
    size_t n_freqs;

    double *acc_re;
    double *acc_im;
} probe_dft_t;

int  probe_dft_create(probe_dft_t *probe, size_t i, size_t j, size_t k, field_component_t component, const double *freqs_hz, size_t n_freqs);

void probe_dft_destroy(probe_dft_t *probe);

void probe_dft_record(probe_dft_t *probe, const fdtd_grid_t *g, double t, double dt);

void probe_dft_result(const probe_dft_t *probe, size_t f_idx, double *out_magnitude, double *out_phase_rad);

#endif