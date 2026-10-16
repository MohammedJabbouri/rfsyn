/*
 * validate_long_term_stability.c -- late-time stability screening.
 *
 * Cases (all 3-D, library update routines):
 *   vacuum_cpml        vacuum, CPML on all faces
 *   eps100_cpml        homogeneous eps_r = 100 (through the CPML)
 *   interface25_cpml   air | eps_r = 25 half-space, interface through CPML
 *   cavity_vacuum      closed PEC box, no CPML  (bounded-energy control)
 *   cavity_interface25 closed PEC box half-filled with eps_r = 25 (control)
 *
 * A differentiated-Gaussian soft source (zero DC) drives Ez and Ex at two
 * off-centre points, then switches off. Every --sample-every steps the
 * program records:
 *   step,t,source_off,energy_diagnostic,max_e,max_h,pml_max_e,pml_max_h
 *
 * energy_diagnostic = 1/2 sum(eps E^2 + mu0 H^2) dx^3 using E and H at
 * their staggered half-step times. It is a diagnostic, not the exactly
 * conserved discrete energy, so it oscillates slightly in closed cavities.
 *
 * Checks
 *   all cases : every sample and the final fields are finite
 *   open cases: final/peak energy <= --resid-tol (default 1e-6), and the
 *               growth screen: after source-off the energy may never rise
 *               above --growth-factor (default 10) x its running minimum
 *               (ignoring values below 1e-20 x peak, i.e. round-off).
 *               Recording continues after a growth flag.
 *   cavities  : |mean(last 10%) / mean(first 10%) - 1| <= --drift-tol
 *               (default 1e-2) over the post-source record.
 * A growth-screen flag requires investigation; it is not by itself proof
 * of exponential instability.
 *
 * Options:
 *   --steps N          (default 100000)
 *   --n N              cells per side (default 24)
 *   --npml N           (default 6)
 *   --case NAME        run only one case
 *   --sample-every N   (default 100)
 *   --strict           nonzero exit if any check fails
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grid.h"
#include "cpml.h"
#include "materials.h"
#include "accuracy_test.h"

typedef struct {
    const char *name;
    int open;            /* 1 = CPML, 0 = PEC cavity */
    double eps_r;        /* value of the filled region */
    int half_space;      /* 1 = only i >= nx/2 filled */
} case_t;

static const case_t CASES[] = {
    {"vacuum_cpml",        1, 1.0,   0},
    {"eps100_cpml",        1, 100.0, 0},
    {"interface25_cpml",   1, 25.0,  1},
    {"cavity_vacuum",      0, 1.0,   0},
    {"cavity_interface25", 0, 25.0,  1},
};
#define N_CASES (sizeof CASES / sizeof CASES[0])

static int in_pml(const fdtd_grid_t *g, size_t i, size_t j, size_t k, size_t npml) {
    return i < npml || j < npml || k < npml ||
           i >= g->nx - npml || j >= g->ny - npml || k >= g->nz - npml;
}

static int run_case(const case_t *c, size_t n, size_t npml, size_t steps, size_t every,
                    double resid_tol, double growth_factor, double drift_tol, FILE *fsum)
{
    const double dx = 1e-3;
    fdtd_grid_t *g = fdtd_grid_create(n, n, n, dx, 0.99);
    if (!g) { fprintf(stderr, "grid alloc failed\n"); return -1; }
    cpml_t *pml = NULL;
    materials_t *mat = NULL;
    if (c->open) {
        pml = cpml_create(g, cpml_default_params(npml));
        if (!pml) { fprintf(stderr, "cpml_create failed\n"); fdtd_grid_destroy(g); return -1; }
    }
    size_t i_fill = c->half_space ? n / 2 : 0;
    if (c->eps_r != 1.0) {
        mat = materials_create(g);
        if (!mat || materials_set_box(mat, g, i_fill, n + 1, 0, n + 1, 0, n + 1, c->eps_r, 0.0) != 0) {
            fprintf(stderr, "materials setup failed\n");
            materials_destroy(mat); cpml_destroy(pml); fdtd_grid_destroy(g); return -1;
        }
    }

    const double dt = g->dt, tau = 8.0 * dt, t0 = 4.5 * tau;
    const size_t src_off_step = (size_t)ceil(2.0 * t0 / dt) + 1;
    size_t s1 = fdtd_index(g, n / 2 - 2, n / 2 + 1, n / 2 + 3);
    size_t s2 = fdtd_index(g, n / 2 + 3, n / 2 - 2, n / 2 - 1);

    char path[96];
    snprintf(path, sizeof path, "stability_%s.csv", c->name);
    FILE *f = at_xfopen(path);
    fprintf(f, "step,t,source_off,energy_diagnostic,max_e,max_h,pml_max_e,pml_max_h\n");

    const size_t total = (g->nx + 1) * (g->ny + 1) * (g->nz + 1);
    size_t n_samp_cap = steps / every + 2, n_samp = 0;
    double *esamp = at_xcalloc(n_samp_cap, sizeof(double), "energy samples");
    size_t *ssamp = at_xcalloc(n_samp_cap, sizeof(size_t), "step samples");

    double peak = 0.0, run_min = INFINITY;
    long growth_step = -1;
    int finite_ok = 1;
    const double cell_vol = dx * dx * dx;

    for (size_t s = 0; s < steps; s++) {
        if (c->open) {
            fdtd_update_h_cpml(g, pml);
            if (mat) fdtd_update_e_cpml_materials(g, pml, mat); else fdtd_update_e_cpml(g, pml);
        } else {
            fdtd_update_h(g);
            if (mat) fdtd_update_e_materials(g, mat); else fdtd_update_e(g);
        }
        double t = (double)(s + 1) * dt;
        if (s < src_off_step) {
            double w = at_gauss_deriv(t, t0, tau);
            g->ez[s1] += w;
            g->ex[s2] += 0.7 * w;
        }

        if (s % every == 0 || s == steps - 1) {
            double en = 0, me = 0, mh = 0, pe = 0, ph = 0;
            for (size_t i = 0; i <= g->nx; i++) {
                double eps = AT_EPS0 * ((mat && i >= i_fill) ? c->eps_r : 1.0);
                for (size_t j = 0; j <= g->ny; j++)
                    for (size_t k = 0; k <= g->nz; k++) {
                        size_t q = fdtd_index(g, i, j, k);
                        double e2 = g->ex[q] * g->ex[q] + g->ey[q] * g->ey[q] + g->ez[q] * g->ez[q];
                        double h2 = g->hx[q] * g->hx[q] + g->hy[q] * g->hy[q] + g->hz[q] * g->hz[q];
                        en += 0.5 * (eps * e2 + AT_MU0 * h2) * cell_vol;
                        double ae = sqrt(e2), ah = sqrt(h2);
                        if (ae > me) me = ae;
                        if (ah > mh) mh = ah;
                        if (c->open && in_pml(g, i, j, k, npml)) {
                            if (ae > pe) pe = ae;
                            if (ah > ph) ph = ah;
                        }
                        if (!isfinite(e2) || !isfinite(h2)) finite_ok = 0;
                    }
            }
            int off = s >= src_off_step;
            fprintf(f, "%zu,%.9e,%d,%.9e,%.9e,%.9e,%.9e,%.9e\n", s, t, off, en, me, mh, pe, ph);
            if (!isfinite(en)) finite_ok = 0;
            if (en > peak) peak = en;
            if (off && n_samp < n_samp_cap) { esamp[n_samp] = en; ssamp[n_samp] = s; n_samp++; }
            if (off && c->open) {
                if (en < run_min) run_min = en;
                if (growth_step < 0 && en > 1e-20 * peak && en > growth_factor * run_min)
                    growth_step = (long)s;   /* flag but keep recording */
            }
            if (!finite_ok) { fprintf(stderr, "%s: non-finite values at step %zu -- stopping\n", c->name, s); break; }
        }
    }
    finite_ok &= at_all_finite(g->ex, total) && at_all_finite(g->ey, total) && at_all_finite(g->ez, total) &&
                 at_all_finite(g->hx, total) && at_all_finite(g->hy, total) && at_all_finite(g->hz, total);
    int io_err = at_xfclose(f, path);

    double final_e = n_samp ? esamp[n_samp - 1] : NAN;
    double metric = NAN;
    int pass = finite_ok && n_samp >= 10 && !io_err;
    const char *metric_name;
    if (c->open) {
        metric_name = "final/peak";
        metric = final_e / peak;
        pass = pass && metric <= resid_tol && growth_step < 0;
    } else {
        metric_name = "drift";
        size_t w = n_samp / 10; if (w < 1) w = 1;
        double a = 0, b = 0;
        for (size_t q = 0; q < w; q++) { a += esamp[q]; b += esamp[n_samp - 1 - q]; }
        metric = (n_samp >= 2 && a > 0) ? fabs(b / a - 1.0) : NAN;
        pass = pass && isfinite(metric) && metric <= drift_tol;
    }

    printf("%-20s %8zu %12.4e %12.4e %-10s %12.4e %10s %s\n", c->name, steps, peak, final_e,
           metric_name, metric, growth_step < 0 ? "none" : "FLAGGED", pass ? "PASS" : "FAIL");
    if (growth_step >= 0) printf("    growth screen first flagged at step %ld\n", growth_step);
    fprintf(fsum, "%s,%d,%g,%zu,%zu,%zu,%.9e,%.9e,%s,%.9e,%ld,%d,%d\n", c->name, c->open, c->eps_r,
            n, npml, steps, peak, final_e, metric_name, metric, growth_step, finite_ok, pass);

    free(esamp); free(ssamp);
    materials_destroy(mat); cpml_destroy(pml); fdtd_grid_destroy(g);
    return pass ? 0 : 1;
}

int main(int argc, char **argv) {
    const int strict = at_has_flag(argc, argv, "--strict");
    const size_t steps = (size_t)at_get_long(argc, argv, "--steps", 100000);
    const size_t n = (size_t)at_get_long(argc, argv, "--n", 24);
    const size_t npml = (size_t)at_get_long(argc, argv, "--npml", 6);
    const size_t every = (size_t)at_get_long(argc, argv, "--sample-every", 100);
    const double resid_tol = at_get_double(argc, argv, "--resid-tol", 1e-6);
    const double growth = at_get_double(argc, argv, "--growth-factor", 10.0);
    const double drift = at_get_double(argc, argv, "--drift-tol", 1e-2);
    const char *only = at_get_str(argc, argv, "--case", NULL);

    if (steps < 1000 || every < 1 || n < 2 * npml + 8 || npml < 2) {
        fprintf(stderr, "bad options: need steps>=1000, n>=2*npml+8, npml>=2\n"); return 2;
    }

    printf("=== validate_long_term_stability ===\n");
    printf("grid %zu^3, npml=%zu, steps=%zu, sample every %zu\n", n, npml, steps, every);
    printf("%-20s %8s %12s %12s %-10s %12s %10s %s\n",
           "case", "steps", "peak_energy", "final_energy", "metric", "value", "growth", "result");

    const char *sum_path = "stability_summary.csv";
    FILE *fs = at_xfopen(sum_path);
    fprintf(fs, "case,open,eps_r,n,npml,steps,peak_energy,final_energy,metric,metric_value,growth_flag_step,finite,pass\n");

    int n_fail = 0, n_err = 0, ran = 0;
    for (size_t c = 0; c < N_CASES; c++) {
        if (only && strcmp(only, CASES[c].name) != 0) continue;
        ran++;
        int r = run_case(&CASES[c], n, npml, steps, every, resid_tol, growth, drift, fs);
        if (r < 0) n_err++; else if (r > 0) n_fail++;
    }
    if (at_xfclose(fs, sum_path)) n_err++;
    if (!ran) { fprintf(stderr, "unknown --case '%s'\n", only); return 2; }
    if (n_err) { printf("\nRESULT: ERROR (execution failure)\n"); return 1; }
    printf("\nRESULT: %d case(s) failed%s\n", n_fail, strict ? "" : " (diagnostic mode; use --strict)");
    return (strict && n_fail) ? 1 : 0;
}
