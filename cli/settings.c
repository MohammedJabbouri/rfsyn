#include "settings.h"
#include "presets.h"
#include "fdtd_stage.h"
#include "../engine/noise/awgn.h"

static void build_fdtd_config(const config_t *cfg, fdtd_sim_config_t *s) {
    s->nx = (size_t)config_get_long(cfg, "fdtd", "nx", PRESET_FDTD_NX);
    s->ny = (size_t)config_get_long(cfg, "fdtd", "ny", PRESET_FDTD_NY);
    s->nz = (size_t)config_get_long(cfg, "fdtd", "nz", PRESET_FDTD_NZ);
    s->dx = config_get_double(cfg, "fdtd", "dx", PRESET_FDTD_DX);
    s->courant_safety = config_get_double(cfg, "fdtd", "courant", PRESET_FDTD_COURANT);
    s->n_steps = (size_t)config_get_long(cfg, "fdtd", "n_steps", PRESET_FDTD_NSTEPS);
    s->npml = (size_t)config_get_long(cfg, "fdtd", "npml", PRESET_FDTD_NPML);
    s->t0_over_dt = config_get_double(cfg, "fdtd", "t0", PRESET_FDTD_T0);
    s->tau_over_dt = config_get_double(cfg, "fdtd", "tau", PRESET_FDTD_TAU);

    s->src_component = FIELD_EZ;
    s->probe_component = FIELD_EZ;

    long si = config_get_long(cfg, "fdtd", "src_i", PRESET_FDTD_SRC_I);
    long sj = config_get_long(cfg, "fdtd", "src_j", PRESET_FDTD_SRC_J);
    long sk = config_get_long(cfg, "fdtd", "src_k", PRESET_FDTD_SRC_K);
    s->src_i = (si < 0) ? s->nx / 2 : (size_t)si;
    s->src_j = (sj < 0) ? s->ny / 2 : (size_t)sj;
    s->src_k = (sk < 0) ? s->nz / 2 : (size_t)sk;

    long pi = config_get_long(cfg, "fdtd", "probe_i", PRESET_FDTD_PROBE_I);
    long pj = config_get_long(cfg, "fdtd", "probe_j", PRESET_FDTD_PROBE_J);
    long pk = config_get_long(cfg, "fdtd", "probe_k", PRESET_FDTD_PROBE_K);
    s->probe_i = (pi < 0) ? s->nx / 2 + 6 : (size_t)pi;
    s->probe_j = (pj < 0) ? s->ny / 2 : (size_t)pj;
    s->probe_k = (pk < 0) ? s->nz / 2 : (size_t)pk;
}

int settings_build_chain(
    const config_t *cfg,
    uint64_t example_index,
    transform_t *stages_out,
    size_t max_stages,
    size_t *n_stages_out
) {
    size_t n = 0;

    if (!n_stages_out) return -1;
    *n_stages_out = 0;

    if (!cfg || !stages_out) return -1;

    if (config_get_bool(cfg, "fdtd", "enabled", PRESET_FDTD_ENABLED)) {
        if (n >= max_stages) return -1;

        fdtd_sim_config_t sim = {0};
        build_fdtd_config(cfg, &sim);

        transform_t stage = fdtd_stage_create(&sim);
        if (!stage.apply) return -1;

        stages_out[n++] = stage;
        *n_stages_out = n;
    }

    if (config_get_bool(cfg, "awgn", "enabled", PRESET_REALISTIC_AWGN_ENABLED)) {
        if (n >= max_stages) return -1;

        double snr_db = config_get_double(cfg, "awgn", "snr_db", PRESET_REALISTIC_AWGN_SNRDB);

        long seed = config_get_long(cfg, "awgn", "seed", PRESET_REALISTIC_AWGN_SEED);

        transform_t stage = awgn_create((float)snr_db, (uint64_t)seed, example_index);

        if (!stage.apply) return -1;

        stages_out[n++] = stage;
        *n_stages_out = n;
    }

    return 0;
}