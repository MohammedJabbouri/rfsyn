#ifndef ENGINE_CORE_RF_SIGNAL_H
#define ENGINE_CORE_RF_SIGNAL_H

#include <complex.h>
#include <stddef.h>

typedef void *(*metadata_copy_fn)(const void *metadata);
typedef void  (*metadata_free_fn)(void *metadata);

typedef struct signal {
	float complex *samples;
	size_t n_samples;
	double sample_rate_hz;
	double center_freq_hz;

	void *metadata;
	metadata_copy_fn meta_copy;
	metadata_free_fn meta_free;
} signal_t;

signal_t *signal_create(size_t n_samples, double sample_rate_hz, double center_freq_hz);

void signal_destroy(signal_t *sig);

int signal_resize(signal_t *sig, size_t new_n_samples);

signal_t *signal_copy(const signal_t *sig);

#endif