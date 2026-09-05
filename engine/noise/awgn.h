#ifndef ENGINE_NOISE_AWGN_H
#define ENGINE_NOISE_AWGN_H

#include <stdint.h>
#include "../core/transform.h"

transform_t awgn_create(float snr_db, uint64_t seed, uint64_t stream);

#endif