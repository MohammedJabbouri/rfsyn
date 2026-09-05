#include "settings.h"
#include "presets.h"
#include "../engine/noise/awgn.h"

int settings_build_chain(const config_t *cfg, uint64_t example_index, transform_t *stages_out, size_t max_stages, size_t *n_stages_out) {
    size_t n = 0;
    if (n_stages_out) *n_stages_out = 0;

    if (config_get_bool(cfg, "awgn", "enabled", PRESET_REALISTIC_AWGN_ENABLED)) {
        if (n >= max_stages) return -1;
        double snr_db = config_get_double(cfg, "awgn", "snr_db", PRESET_REALISTIC_AWGN_SNRDB);
        long seed = config_get_long(cfg, "awgn", "seed", PRESET_REALISTIC_AWGN_SEED);
        stages_out[n++] = awgn_create((float)snr_db, (uint64_t)seed, example_index);
        if (n_stages_out) *n_stages_out = n;
    }

    return 0;
}