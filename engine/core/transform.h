#ifndef ENGINE_CORE_TRANSFORM_H
#define ENGINE_CORE_TRANSFORM_H

#include "rf_signal.h"

typedef struct {
    int (*apply)(void *ctx, signal_t *sig);
    void (*destroy)(void *ctx);
    void *ctx;
} transform_t;

#endif