#include "chain.h"

int chain_apply(const transform_t *stages, size_t n_stages, signal_t *sig) {
    if (!stages || !sig) return -1;

    for (size_t i = 0; i < n_stages; i++) {
        if (!stages[i].apply) continue;
        if (stages[i].apply(stages[i].ctx, sig) != 0) {
            return -1;
        }
    }
    return 0;
}

void chain_free_stages(transform_t *stages, size_t n_stages) {
    if (!stages) return;
    for (size_t i = 0; i < n_stages; i++) {
        if (stages[i].destroy) {
            stages[i].destroy(stages[i].ctx);
        }
    }
}

// WILL REPLACE THIS LATER ON