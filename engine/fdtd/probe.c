#include "probe.h"
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const double *field_ptr(const fdtd_grid_t *g, field_component_t c) {
    switch (c) {
        case FIELD_EX: return g->ex;
        case FIELD_EY: return g->ey;
        case FIELD_EZ: return g->ez;
        case FIELD_HX: return g->hx;
        case FIELD_HY: return g->hy;
        case FIELD_HZ: return g->hz;
        default: return NULL;
    }
}

int probe_td_create(probe_td_t *probe, size_t i, size_t j, size_t k,
                     field_component_t component, size_t n_steps)
{
    if (!probe || n_steps == 0) return -1;
    probe->buffer = malloc(n_steps * sizeof(double));
    if (!probe->buffer) return -1;
    probe->i = i; probe->j = j; probe->k = k;
    probe->component = component;
    probe->n_steps = n_steps;
    probe->next_index = 0;
    return 0;
}

void probe_td_destroy(probe_td_t *probe) {
    if (!probe) return;
    free(probe->buffer);
    probe->buffer = NULL;
}

void probe_td_record(probe_td_t *probe, const fdtd_grid_t *g) {
    if (!probe || !g || !probe->buffer) return;
    if (probe->next_index >= probe->n_steps) return;

    const double *field = field_ptr(g, probe->component);
    if (!field) return;

    size_t idx = fdtd_index(g, probe->i, probe->j, probe->k);
    probe->buffer[probe->next_index++] = field[idx];
}

int probe_dft_create(probe_dft_t *probe, size_t i, size_t j, size_t k, field_component_t component, const double *freqs_hz, size_t n_freqs)
{
    if (!probe || !freqs_hz || n_freqs == 0) return -1;

    probe->acc_re = calloc(n_freqs, sizeof(double));
    probe->acc_im = calloc(n_freqs, sizeof(double));
    if (!probe->acc_re || !probe->acc_im) {
        free(probe->acc_re);
        free(probe->acc_im);
        return -1;
    }

    probe->i = i; probe->j = j; probe->k = k;
    probe->component = component;
    probe->freqs_hz = freqs_hz;
    probe->n_freqs = n_freqs;
    return 0;
}

void probe_dft_destroy(probe_dft_t *probe) {
    if (!probe) return;
    free(probe->acc_re);
    free(probe->acc_im);
    probe->acc_re = NULL;
    probe->acc_im = NULL;
}

void probe_dft_record(probe_dft_t *probe, const fdtd_grid_t *g, double t, double dt) {
    if (!probe || !g) return;

    const double *field = field_ptr(g, probe->component);
    if (!field) return;

    size_t idx = fdtd_index(g, probe->i, probe->j, probe->k);
    double val = field[idx];

    for (size_t f = 0; f < probe->n_freqs; f++) {
        double omega = 2.0 * M_PI * probe->freqs_hz[f];
        probe->acc_re[f] += val * cos(omega * t) * dt;
        probe->acc_im[f] -= val * sin(omega * t) * dt;
    }
}

void probe_dft_result(const probe_dft_t *probe, size_t f_idx, double *out_magnitude, double *out_phase_rad)
{
    if (!probe || f_idx >= probe->n_freqs) {
        if (out_magnitude) *out_magnitude = 0.0;
        if (out_phase_rad) *out_phase_rad = 0.0;
        return;
    }
    double re = probe->acc_re[f_idx];
    double im = probe->acc_im[f_idx];
    if (out_magnitude) *out_magnitude = sqrt(re * re + im * im);
    if (out_phase_rad) *out_phase_rad = atan2(im, re);
}