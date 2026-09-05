#ifndef ENGINE_CORE_TRANSFORM_H
#define ENGINE_CORE_TRANSFORM_H

#include "signal.h"

typedef struct {
    int (*apply)(void *ctx, signal_t *sig);
    void (*destroy)(void *ctx);
    void *ctx;
} transform_t;

#endif