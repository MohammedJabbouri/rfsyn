#include "rf_signal.h"
#include <stdlib.h>
#include <string.h>
#include <complex.h>

signal_t *signal_create(size_t n_samples, double sample_rate_hz, double center_freq_hz) {
	signal_t *sig = malloc(sizeof(signal_t));
	if (!sig) return NULL;

	sig->samples = calloc(n_samples,sizeof(float complex));
	if (!sig->samples && n_samples > 0) {
		free(sig);
		return NULL;
	}
	sig->n_samples = n_samples;
	sig->sample_rate_hz = sample_rate_hz;
	sig->center_freq_hz = center_freq_hz;
	sig->metadata = NULL;
	sig->meta_copy = NULL;
	sig->meta_free = NULL;

	return sig;
}

void signal_destroy(signal_t * sig) {
	if (!sig) return;
	free(sig->samples);
	free(sig);
}

int signal_resize(signal_t *sig, size_t new_n_samples) {
	if (!sig) return -1;

	float complex *grown = realloc(sig->samples, new_n_samples * sizeof(float complex));
	if (!grown && new_n_samples > 0) {
		return -1;
	}

	size_t old_n_samples = sig->n_samples;
	sig->samples   = grown;
	sig->n_samples = new_n_samples;

	if (new_n_samples > old_n_samples) {
		size_t old_bytes = old_n_samples * sizeof(float complex);
		size_t new_bytes = new_n_samples * sizeof(float complex);
		memset((char *)grown + old_bytes, 0, new_bytes - old_bytes);
	}

	return 0;
}

signal_t *signal_copy(const signal_t *sig) {
	if (!sig) return NULL;

	signal_t *copy = signal_create(sig->n_samples, sig->sample_rate_hz, sig->center_freq_hz);
	if (!copy) return NULL;

	memcpy(copy->samples, sig->samples, sig->n_samples * sizeof(float complex));

	copy->metadata = sig->metadata;
	
	return copy;
}