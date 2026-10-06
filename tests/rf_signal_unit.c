#include <stdio.h>
#include <string.h>
#include <complex.h>
#include "../engine/core/rf_signal.h"

static int check(int cond, const char *what) {
    if (!cond) {
        fprintf(stderr, "FAIL: %s\n", what);
    }
    return cond;
}

int main(void) {
    int ok = 1;

    signal_t *sig = signal_create(100, 1e6, 915e6);
    ok &= check(sig != NULL, "signal_create returned non-NULL");
    ok &= check(sig->n_samples == 100, "initial n_samples == 100");
    ok &= check(sig->meta_copy == NULL, "meta_copy initialized to NULL");
    ok &= check(sig->meta_free == NULL, "meta_free initialized to NULL");

    for (size_t i = 0; i < sig->n_samples; i++) {
        sig->samples[i] = (float)i + 0.0f * I;
    }

    ok &= check(signal_resize(sig, 200) == 0, "resize grow 100 -> 200 succeeds");
    ok &= check(sig->n_samples == 200, "n_samples updated after grow");
    for (size_t i = 0; i < 100; i++) {
        ok &= check(crealf(sig->samples[i]) == (float)i, "grow preserves prefix value");
    }
    for (size_t i = 100; i < 200; i++) {
        ok &= check(sig->samples[i] == 0.0f, "grow zero-fills new tail");
    }

    ok &= check(signal_resize(sig, 50) == 0, "resize shrink 200 -> 50 succeeds");
    ok &= check(sig->n_samples == 50, "n_samples updated after shrink");
    for (size_t i = 0; i < 50; i++) {
        sig->samples[i] = 42.0f + 0.0f * I;
    }
    for (size_t i = 0; i < 50; i++) {
        ok &= check(sig->samples[i] == 42.0f, "post-shrink write/read round-trips");
    }

    ok &= check(signal_resize(sig, 50) == 0, "resize no-op 50 -> 50 succeeds");
    ok &= check(sig->n_samples == 50, "n_samples unchanged after no-op resize");
    sig->samples[49] = 7.0f + 0.0f * I;
    ok &= check(sig->samples[49] == 7.0f, "post-no-op-resize write/read round-trips");

    ok &= check(signal_resize(sig, 0) == 0, "resize to 0 succeeds");
    ok &= check(sig->n_samples == 0, "n_samples == 0 after resize to 0");

    signal_destroy(sig);

    signal_t *a = signal_create(10, 1e6, 915e6);
    for (size_t i = 0; i < 10; i++) a->samples[i] = (float)i + 1.0f * I;
    signal_t *b = signal_copy(a);
    ok &= check(b != NULL, "signal_copy returned non-NULL");
    ok &= check(b->samples != a->samples, "copy has its own sample buffer");
    for (size_t i = 0; i < 10; i++) {
        ok &= check(b->samples[i] == a->samples[i], "copy has matching values");
    }
    signal_destroy(a);

    ok &= check(b->samples[3] == 3.0f + 1.0f * I, "copy survives original's destruction");
    signal_destroy(b);

    if (!ok) {
        fprintf(stderr, "\nFAIL\n");
        return 1;
    }
    printf("PASS\n");
    return 0;
}