#include "fdtd_stage.h"
#include "../engine/core/rf_signal.h"
#include <stdlib.h>
#include <complex.h>

typedef struct {
    fdtd_sim_config_t cfg;
} fdtd_stage_ctx_t;

static int fdtd_stage_apply(void *ctx_v, signal_t *sig) {
    fdtd_stage_ctx_t *ctx = (fdtd_stage_ctx_t *)ctx_v;
    if (!ctx || !sig) return -1;

    if (signal_resize(sig, ctx->cfg.n_steps) != 0) return -1;

    double *buf = malloc(ctx->cfg.n_steps * sizeof(double));
    if (!buf) return -1;

    double fs = 0.0;
    int rc = fdtd_simulate_echo(&ctx->cfg, buf, &fs);
    if (rc != 0) { free(buf); return -1; }

    for (size_t i = 0; i < sig->n_samples; i++)
        sig->samples[i] = (float)buf[i] + 0.0f * I;

    free(buf);
    sig->sample_rate_hz = fs;
    return 0;
}

static void fdtd_stage_destroy(void *ctx_v) {
    free(ctx_v);
}

transform_t fdtd_stage_create(const fdtd_sim_config_t *cfg) {
    if (!cfg) return (transform_t){ NULL, NULL, NULL };
    fdtd_stage_ctx_t *ctx = malloc(sizeof(fdtd_stage_ctx_t));
    if (!ctx) return (transform_t){ NULL, NULL, NULL };
    ctx->cfg = *cfg;
    return (transform_t){ fdtd_stage_apply, fdtd_stage_destroy, ctx };
}