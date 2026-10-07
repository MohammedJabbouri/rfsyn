#include "presets.h"
#include <string.h>

config_t *preset_build(const char *name) {
    if (!name || strcmp(name, "realistic") != 0) return NULL;

    config_t *cfg = config_create_empty();
    if (!cfg) return NULL;

    config_set_long(cfg, "job", "count", PRESET_REALISTIC_JOB_COUNT);
    config_set_long(cfg, "job", "n_samples", PRESET_REALISTIC_JOB_NSAMPLES);
    config_set_string(cfg, "job", "output_dir", PRESET_REALISTIC_JOB_OUTPUTDIR);

    config_set_bool(cfg, "awgn", "enabled", PRESET_REALISTIC_AWGN_ENABLED);
    config_set_double(cfg, "awgn", "snr_db", PRESET_REALISTIC_AWGN_SNRDB);
    config_set_long(cfg, "awgn", "seed", PRESET_REALISTIC_AWGN_SEED);

    config_set_bool(cfg, "fdtd", "enabled", PRESET_FDTD_ENABLED);
    config_set_long(cfg, "fdtd", "nx", PRESET_FDTD_NX);
    config_set_long(cfg, "fdtd", "ny", PRESET_FDTD_NY);
    config_set_long(cfg, "fdtd", "nz", PRESET_FDTD_NZ);
    config_set_double(cfg, "fdtd", "dx", PRESET_FDTD_DX);
    config_set_double(cfg, "fdtd", "courant", PRESET_FDTD_COURANT);
    config_set_long(cfg, "fdtd", "n_steps", PRESET_FDTD_NSTEPS);
    config_set_long(cfg, "fdtd", "npml", PRESET_FDTD_NPML);
    config_set_long(cfg, "fdtd", "src_i", PRESET_FDTD_SRC_I);
    config_set_long(cfg, "fdtd", "src_j", PRESET_FDTD_SRC_J);
    config_set_long(cfg, "fdtd", "src_k", PRESET_FDTD_SRC_K);
    config_set_long(cfg, "fdtd", "probe_i", PRESET_FDTD_PROBE_I);
    config_set_long(cfg, "fdtd", "probe_j", PRESET_FDTD_PROBE_J);
    config_set_long(cfg, "fdtd", "probe_k", PRESET_FDTD_PROBE_K);
    config_set_double(cfg, "fdtd", "t0", PRESET_FDTD_T0);
    config_set_double(cfg, "fdtd", "tau", PRESET_FDTD_TAU);

    config_set_bool(cfg, "material_catalog", "enabled", 0);
    config_set_double(cfg, "material_catalog", "evaluation_frequency_hz", 1e9);
    config_set_bool(cfg, "material_catalog", "allow_proxy", 0);

    return cfg;
}