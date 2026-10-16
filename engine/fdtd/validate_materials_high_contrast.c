#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "accuracy_test.h"

typedef struct {
    double R, T, R_th, T_th;
    double cells_per_wl_high;
    double dx;
    size_t n_cells, n_steps;
} hc_result_t;

static int run_case(double eps1, double eps2, int scale, double f0, int staircase,
                    const char *wave_path, hc_result_t *out)
{
    const double dx = 1e-3 / (double)scale;
    const double dt = 0.99 * dx / (AT_C0 * sqrt(3.0));
    const double n1 = sqrt(eps1), n2 = sqrt(eps2);
    const double v1 = AT_C0 / n1, v2 = AT_C0 / n2;

    const double tau = 1.0 / f0;
    const double t0 = 4.5 * tau;

    const double d_src = 0.120, d1 = 0.060, d2 = 0.060;
    double t_end = 2.0 * t0 + fmax((d_src + d1) / v1, d_src / v1 + d2 / v2) + 2.0 * tau;
    size_t n_steps = (size_t)ceil(t_end / dt);

    size_t pad = (size_t)ceil(0.55 * AT_C0 * t_end / dx) + 20;
    size_t c_src = (size_t)llround(d_src / dx), c1 = (size_t)llround(d1 / dx), c2 = (size_t)llround(d2 / dx);
    size_t is = pad;
    size_t ii = is + c_src;
    size_t p1 = ii - c1, p2 = ii + c2;
    size_t N = p2 + pad;

    double *ea = at_xcalloc(N + 1, sizeof(double), "ea"), *ha = at_xcalloc(N, sizeof(double), "ha");
    double *eb = at_xcalloc(N + 1, sizeof(double), "eb"), *hb = at_xcalloc(N, sizeof(double), "hb");
    double *cba = at_xcalloc(N + 1, sizeof(double), "cba"), *cbb = at_xcalloc(N + 1, sizeof(double), "cbb");
    const double chh = dt / (AT_MU0 * dx);
    for (size_t i = 0; i <= N; i++) {
        double eb_r = (i < ii) ? eps1 : (i > ii) ? eps2 : (staircase ? eps2 : 0.5 * (eps1 + eps2));
        cba[i] = dt / (AT_EPS0 * eps1 * dx);
        cbb[i] = dt / (AT_EPS0 * eb_r * dx);
    }

    FILE *fw = wave_path ? at_xfopen(wave_path) : NULL;
    if (fw) fprintf(fw, "step,t,incident_p1,reflected_p1,transmitted_p2\n");
    at_dft_t Di = {0, 0}, Dr = {0, 0}, Dt = {0, 0};
    int rc = 0;

    for (size_t n = 0; n < n_steps; n++) {
        for (size_t i = 0; i < N; i++) {
            ha[i] += chh * (ea[i + 1] - ea[i]);
            hb[i] += chh * (eb[i + 1] - eb[i]);
        }
        for (size_t i = 1; i < N; i++) {
            ea[i] += cba[i] * (ha[i] - ha[i - 1]);
            eb[i] += cbb[i] * (hb[i] - hb[i - 1]);
        }
        double t = (double)(n + 1) * dt;
        double s = at_mod_gauss(t, f0, t0, tau);
        ea[is] += s; eb[is] += s;

        double inc = ea[p1], refl = eb[p1] - ea[p1], tr = eb[p2];
        if (!isfinite(inc) || !isfinite(refl) || !isfinite(tr)) {
            fprintf(stderr, "non-finite field at step %zu\n", n); rc = 1; break;
        }
        at_dft_add(&Di, inc, t, f0, dt);
        at_dft_add(&Dr, refl, t, f0, dt);
        at_dft_add(&Dt, tr, t, f0, dt);
        if (fw && (n % (size_t)scale == 0))
            fprintf(fw, "%zu,%.9e,%.9e,%.9e,%.9e\n", n, t, inc, refl, tr);
    }
    if (fw && at_xfclose(fw, wave_path)) rc = 1;

    double gi = at_dft_abs2(Di);
    out->R = at_dft_abs2(Dr) / gi;
    out->T = (n2 / n1) * at_dft_abs2(Dt) / gi;
    double r = (n1 - n2) / (n1 + n2);
    out->R_th = r * r;
    out->T_th = 1.0 - r * r;
    out->dx = dx;
    out->cells_per_wl_high = (AT_C0 / f0) / (fmax(n1, n2) * dx);
    out->n_cells = N;
    out->n_steps = n_steps;

    free(ea); free(ha); free(eb); free(hb); free(cba); free(cbb);
    return rc;
}

int main(int argc, char **argv) {
    const int strict = at_has_flag(argc, argv, "--strict");
    const int staircase = at_has_flag(argc, argv, "--staircase");
    const int waves = !at_has_flag(argc, argv, "--no-waveforms");
    const int scale = (int)at_get_long(argc, argv, "--scale", 1);
    const double f0 = at_get_double(argc, argv, "--freq", 1e9);
    const double tol = at_get_double(argc, argv, "--tol", 0.02);
    const char *dir = at_get_str(argc, argv, "--direction", "both");

    if (scale != 1 && scale != 2 && scale != 4) { fprintf(stderr, "--scale must be 1, 2 or 4\n"); return 2; }
    if (f0 <= 0) { fprintf(stderr, "--freq must be positive\n"); return 2; }
    int do_in = !strcmp(dir, "both") || !strcmp(dir, "in");
    int do_out = !strcmp(dir, "both") || !strcmp(dir, "out");
    if (!do_in && !do_out) { fprintf(stderr, "--direction must be in, out or both\n"); return 2; }

    double eps_list[8]; size_t n_eps = 0;
    if (at_get_str(argc, argv, "--eps", NULL)) {
        eps_list[n_eps++] = at_get_double(argc, argv, "--eps", 1.0);
        if (eps_list[0] < 1.0) { fprintf(stderr, "--eps must be >= 1\n"); return 2; }
    } else {
        const double e[] = {1, 4, 10, 25, 50, 100};
        for (size_t i = 0; i < 6; i++) eps_list[n_eps++] = e[i];
    }

    printf("=== validate_materials_high_contrast (1-D Yee kernel) ===\n");
    printf("f0=%.4g GHz scale=%d dx=%.4g mm interface=%s tol=%.3g\n",
           f0 / 1e9, scale, 1.0 / scale, staircase ? "staircase" : "averaged", tol);
    printf("%7s %-14s %8s %10s %10s %10s %10s %10s %9s  %s\n",
           "eps_r", "direction", "cells/wl", "R", "R_theory", "T", "T_theory", "R+T", "|dR|", "result");

    const char *sum_path = "materials_high_contrast_summary.csv";
    FILE *fs = at_xfopen(sum_path);
    fprintf(fs, "eps_r,direction,scale,dx_m,freq_hz,interface,cells_per_wl_high,R,R_theory,T,T_theory,R_plus_T,abs_err_R,abs_err_T,abs_err_closure,pass\n");

    int n_fail = 0, rc = 0;
    for (size_t e = 0; e < n_eps; e++) {
        for (int d = 0; d < 2; d++) {
            if ((d == 0 && !do_in) || (d == 1 && !do_out)) continue;
            double eps1 = d == 0 ? 1.0 : eps_list[e];
            double eps2 = d == 0 ? eps_list[e] : 1.0;
            const char *dname = d == 0 ? "air->dielectric" : "dielectric->air";
            char wpath[96];
            snprintf(wpath, sizeof wpath, "materials_high_contrast_eps%g_%s_s%d.csv",
                     eps_list[e], d == 0 ? "in" : "out", scale);
            hc_result_t r;
            if (run_case(eps1, eps2, scale, f0, staircase, waves ? wpath : NULL, &r)) { rc = 1; continue; }
            double eR = fabs(r.R - r.R_th), eT = fabs(r.T - r.T_th), eC = fabs(r.R + r.T - 1.0);
            int pass = isfinite(r.R) && isfinite(r.T) && eR <= tol && eT <= tol && eC <= tol;
            if (!pass) n_fail++;
            printf("%7g %-14s %8.1f %10.6f %10.6f %10.6f %10.6f %10.6f %9.2e  %s\n",
                   eps_list[e], dname, r.cells_per_wl_high, r.R, r.R_th, r.T, r.T_th, r.R + r.T, eR,
                   pass ? "PASS" : "FAIL");
            fprintf(fs, "%g,%s,%d,%.6e,%.6e,%s,%.4f,%.9e,%.9e,%.9e,%.9e,%.9e,%.6e,%.6e,%.6e,%d\n",
                    eps_list[e], d == 0 ? "in" : "out", scale, r.dx, f0,
                    staircase ? "staircase" : "averaged", r.cells_per_wl_high,
                    r.R, r.R_th, r.T, r.T_th, r.R + r.T, eR, eT, eC, pass);
        }
    }
    if (at_xfclose(fs, sum_path)) rc = 1;

    if (rc) { printf("\nRESULT: ERROR (execution failure)\n"); return 1; }
    printf("\nRESULT: %d case(s) outside tolerance%s\n", n_fail,
           strict ? "" : " (diagnostic mode; use --strict)");
    return (strict && n_fail) ? 1 : 0;
}