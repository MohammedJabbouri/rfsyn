#ifndef ENGINE_CORE_CHAIN_H
#define ENGINE_CORE_CHAIN_H

#include <stddef.h>
#include "signal.h"
#include "transform.h"

int chain_apply(const transform_t *stages, size_t n_stages, signal_t *sig);
void chain_free_stages(transform_t *stages, size_t n_stages);

#endif