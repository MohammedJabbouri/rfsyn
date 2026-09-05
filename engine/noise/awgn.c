#include "awgn.h"
#include "../core/rng.h"
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    float snr_db;
    rng_state_t rng;
} awgn_ctx_t;

static int awgn_apply(void *ctx_v, signal_t *sig) {
    awgn_ctx_t *ctx = (awgn_ctx_t *)ctx_v;
    if (!ctx || !sig || sig->n_samples == 0) return 0;

    double signal_power = 0.0;
    for (size_t i = 0; i < sig->n_samples; i++) {
        float mag = cabsf(sig->samples[i]);
        signal_power += (double)mag * (double)mag;
    }
    signal_power /= (double)sig->n_samples;

    double snr_linear = pow(10.0, (double)ctx->snr_db / 10.0);
    double noise_power = (snr_linear > 0.0) ? signal_power / snr_linear : 0.0;
    double per_axis_std = sqrt(noise_power / 2.0);

    for (size_t i = 0; i < sig->n_samples; i++) {
        float u1 = rng_uniform(&ctx->rng);
        float u2 = rng_uniform(&ctx->rng);
        if (u1 <= 0.0f) u1 = 1e-9f;

        float r = sqrtf(-2.0f * logf(u1));
        float noise_re = (float)per_axis_std * r * cosf(2.0f * (float)M_PI * u2);
        float noise_im = (float)per_axis_std * r * sinf(2.0f * (float)M_PI * u2);

        sig->samples[i] += noise_re + noise_im * I;
    }

    return 0;
}

static void awgn_destroy(void *ctx_v) {
    free(ctx_v);
}

transform_t awgn_create(float snr_db, uint64_t seed, uint64_t stream) {
    awgn_ctx_t *ctx = malloc(sizeof(awgn_ctx_t));
    if (!ctx) {
        return (transform_t){ NULL, NULL, NULL };
    }

    ctx->snr_db = snr_db;
    rng_seed(&ctx->rng, seed, stream);

    return (transform_t){ awgn_apply, awgn_destroy, ctx };
}