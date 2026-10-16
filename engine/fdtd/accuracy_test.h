/*
 * accuracy_test.h -- shared helpers for the new accuracy validators.
 *
 * Header-only (all functions are static) so each validator stays a single
 * translation unit with its own main(). Uses only the public FDTD API
 * (grid.h / cpml.h / materials.h) that is common to the original and the
 * replacement engine files.
 */
#ifndef ENGINE_FDTD_ACCURACY_TEST_H
#define ENGINE_FDTD_ACCURACY_TEST_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if defined(__GNUC__) || defined(__clang__)
#define AT_UNUSED __attribute__((unused))
#else
#define AT_UNUSED
#endif

#define AT_C0   299792458.0
#define AT_EPS0 8.8541878128e-12
#define AT_MU0  1.25663706212e-6

/* ---------------- command-line helpers ---------------- */

static AT_UNUSED int at_has_flag(int argc, char **argv, const char *flag) {
    for (int i = 1; i < argc; i++)
        if (strcmp(argv[i], flag) == 0) return 1;
    return 0;
}

static AT_UNUSED const char *at_get_str(int argc, char **argv, const char *key, const char *def) {
    for (int i = 1; i < argc - 1; i++)
        if (strcmp(argv[i], key) == 0) return argv[i + 1];
    return def;
}

static AT_UNUSED double at_get_double(int argc, char **argv, const char *key, double def) {
    const char *s = at_get_str(argc, argv, key, NULL);
    if (!s) return def;
    char *end = NULL;
    double v = strtod(s, &end);
    if (end == s || *end != '\0') {
        fprintf(stderr, "invalid numeric value for %s: '%s'\n", key, s);
        exit(2);
    }
    return v;
}

static AT_UNUSED long at_get_long(int argc, char **argv, const char *key, long def) {
    const char *s = at_get_str(argc, argv, key, NULL);
    if (!s) return def;
    char *end = NULL;
    long v = strtol(s, &end, 10);
    if (end == s || *end != '\0') {
        fprintf(stderr, "invalid integer value for %s: '%s'\n", key, s);
        exit(2);
    }
    return v;
}

/* ---------------- waveforms ---------------- */

/* Differentiated Gaussian: zero DC content, spectral peak at 1/(sqrt(2)*pi*tau). */
static AT_UNUSED double at_gauss_deriv(double t, double t0, double tau) {
    double a = (t - t0) / tau;
    return -2.0 * a * exp(-a * a);
}

static AT_UNUSED double at_gauss_deriv_peak_freq(double tau) {
    return 1.0 / (sqrt(2.0) * M_PI * tau);
}

/* Gaussian-modulated sinusoid centred on f0. */
static AT_UNUSED double at_mod_gauss(double t, double f0, double t0, double tau) {
    double a = (t - t0) / tau;
    return exp(-a * a) * sin(2.0 * M_PI * f0 * (t - t0));
}

/* ---------------- running single-frequency DFT ---------------- */

typedef struct { double re, im; } at_dft_t;

static AT_UNUSED void at_dft_add(at_dft_t *d, double value, double t, double f, double dt) {
    double w = 2.0 * M_PI * f * t;
    d->re += value * cos(w) * dt;
    d->im -= value * sin(w) * dt;
}

static AT_UNUSED double at_dft_abs2(at_dft_t d) { return d.re * d.re + d.im * d.im; }

/* ---------------- misc ---------------- */

static AT_UNUSED int at_all_finite(const double *a, size_t n) {
    for (size_t i = 0; i < n; i++)
        if (!isfinite(a[i])) return 0;
    return 1;
}

static AT_UNUSED double at_db20(double ratio) {
    return (ratio > 0.0) ? 20.0 * log10(ratio) : -400.0;
}

static AT_UNUSED void *at_xcalloc(size_t n, size_t sz, const char *what) {
    void *p = calloc(n ? n : 1, sz);
    if (!p) {
        fprintf(stderr, "allocation failed: %s (%zu x %zu bytes)\n", what, n, sz);
        exit(1);
    }
    return p;
}

static AT_UNUSED FILE *at_xfopen(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "cannot open '%s' for writing\n", path);
        exit(1);
    }
    return f;
}

static AT_UNUSED int at_xfclose(FILE *f, const char *path) {
    if (ferror(f) || fclose(f) != 0) {
        fprintf(stderr, "write error on '%s'\n", path);
        return 1;
    }
    return 0;
}

/* Approximate bytes needed for one grid + CPML (6 field + 12 psi arrays). */
static AT_UNUSED double at_grid_cpml_mb(size_t nx, size_t ny, size_t nz, int with_cpml) {
    double cells = (double)(nx + 1) * (double)(ny + 1) * (double)(nz + 1);
    return cells * 8.0 * (with_cpml ? 18.0 : 6.0) / (1024.0 * 1024.0);
}

#endif
