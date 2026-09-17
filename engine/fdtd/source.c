#include "source.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double source_waveform_gaussian_derivative(double t, const void *params) {
    const source_gaussian_deriv_params_t *p = (const source_gaussian_deriv_params_t *)params;
    double arg = (t - p->t0) / p->tau;
    return -2.0 * arg * exp(-arg * arg);
}

double source_waveform_ricker(double t, const void *params) {
    const source_ricker_params_t *p = (const source_ricker_params_t *)params;
    double a = M_PI * p->f0_hz * (t - p->t0);
    double a2 = a * a;
    return (1.0 - 2.0 * a2) * exp(-a2);
}

double source_waveform_modulated_gaussian(double t, const void *params) {
    const source_modulated_gaussian_params_t *p = (const source_modulated_gaussian_params_t *)params;
    double arg = (t - p->t0) / p->tau;
    double envelope = exp(-arg * arg);
    return envelope * sin(2.0 * M_PI * p->freq_hz * (t - p->t0));
}

source_t source_make_point(size_t i, size_t j, size_t k, field_component_t component, source_mode_t mode, double amplitude, source_waveform_fn waveform, const void *params)
{
    source_t s;
    s.i = i; s.j = j; s.k = k;
    s.component = component;
    s.mode = mode;
    s.amplitude = amplitude;
    s.waveform = waveform;
    s.params = params;
    return s;
}

void source_apply(const source_t *src, fdtd_grid_t *g, double t) {
    if (!src || !g || !src->waveform) return;
    if (src->i > g->nx || src->j > g->ny || src->k > g->nz) return;

    double *field;
    switch (src->component) {
        case FIELD_EX: field = g->ex; break;
        case FIELD_EY: field = g->ey; break;
        case FIELD_EZ: field = g->ez; break;
        default: return;
    }

    double val = src->amplitude * src->waveform(t, src->params);
    size_t idx = fdtd_index(g, src->i, src->j, src->k);

    if (src->mode == SRC_SOFT) field[idx] += val;
    else                       field[idx] = val;
}