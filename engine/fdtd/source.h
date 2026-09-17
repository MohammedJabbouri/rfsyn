#ifndef ENGINE_FDTD_SOURCE_H
#define ENGINE_FDTD_SOURCE_H

#include <stddef.h>
#include "grid.h"
#include "field_component.h"

typedef enum {
	SRC_SOFT,
	SRC_HARD
} source_mode_t;

typedef double (*source_waveform_fn)(double t, const void *params);

typedef struct {
    size_t i, j, k;
    field_component_t component;  // MUST BE E AS DRIVING H DIRECTLY IS NOT POSSIBLE
    source_mode_t mode;
    double amplitude;
    source_waveform_fn waveform;
    const void *params;
} source_t;

typedef struct { double t0, tau; } source_gaussian_deriv_params_t;
double source_waveform_gaussian_derivative(double t, const void *params);

typedef struct { double f0_hz, t0; } source_ricker_params_t;
double source_waveform_ricker(double t, const void *params);

typedef struct { double freq_hz, t0, tau; } source_modulated_gaussian_params_t;
double source_waveform_modulated_gaussian(double t, const void *params);

source_t source_make_point(size_t i, size_t j, size_t k, field_component_t component, source_mode_t mode, double amplitude, source_waveform_fn waveform, const void *params);

void source_apply(const source_t *src, fdtd_grid_t *g, double t);

#endif