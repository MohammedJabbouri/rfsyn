/*
 * validate_cpml_oblique.c -- CPML reflection versus angle of incidence.
 *
 * Method (image-source geometry):
 *   A soft point dipole sits ds cells in front of the +x CPML face. A row of
 *   probes sits dp cells in front of the same face, offset laterally in y.
 *   The wave reflected by the +x CPML that reaches a probe strikes that face
 *   at   theta = atan(dy / (ds + dp))   from the face normal.
 *   The same simulation is repeated in a reference grid enlarged by
 *   `extra` cells on every side, so (within the recorded window) nothing
 *   reflected by the reference boundary can reach the probes. The
 *   difference test - reference at each probe is the boundary-induced error.
 *
 * Polarizations (incidence plane = x-y):
 *   TE : Ez dipole, Ez probes (E perpendicular to the incidence plane)
 *   TM : Ey dipole, Ey probes (E in the incidence plane)
 *
 * Metrics per probe:
 *   max_norm  = max_t |test-ref| / max_t |ref|        (also in dB)
 *   rel_l2    = ||test-ref||_2 / ||ref||_2
 *   spec_db   = 20 log10 |DFT(test-ref)| / |DFT(ref)| at the pulse peak freq
 *
 * Note: the probe also sees reflections from the other CPML faces (the
 * -x, +-y, +-z faces at larger path lengths). The angle column describes
 * the dominant +x-face reflection only.
 *
 * Options:
 *   --angle DEG        run a single angle instead of the default list
 *   --extended         add 75 and 85 degrees (uses ds=dp=3 unless given)
 *   --npml N           CPML thickness in cells (default 10)
 *   --ds N, --dp N     source / probe standoff from the CPML face (default 6)
 *   --reference-extra N  override the automatic reference enlargement
 *   --pol te|tm|both   (default both)
 *   --gate-db X        strict-mode gate on max_norm (default -30 dB)
 *   --max-mb X         refuse to allocate more than X MB (default 6000)
 *   --strict           return nonzero if any probe exceeds the gate
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "cpml.h"
#include "accuracy_test.h"

#define MAX_ANGLES 16

typedef struct {
    double requested_deg;
    double actual_deg;
    size_t dy;
} angle_t;

typedef enum { POL_TE = 0, POL_TM = 1 } pol_t;

static double *field_of(fdtd_grid_t *g, pol_t pol) {
    return pol == POL_TE ? g->ez : g->ey;
}

/* Runs one simulation and records probe histories into hist[p*n_steps+n]. */
static int run_sim(size_t nx, size_t ny, size_t nz, double dx, size_t npml,
                   size_t off, size_t is, size_t js, size_t ks,
                   size_t ip, const size_t *jp, size_t n_probes,
                   pol_t pol, size_t n_steps, double t0, double tau,
                   double *hist)
{
    fdtd_grid_t *g = fdtd_grid_create(nx, ny, nz, dx, 0.99);
    if (!g) { fprintf(stderr, "fdtd_grid_create failed (%zux%zux%zu)\n", nx, ny, nz); return 1; }
    cpml_t *pml = cpml_create(g, cpml_default_params(npml));
    if (!pml) { fprintf(stderr, "cpml_create failed\n"); fdtd_grid_destroy(g); return 1; }

    size_t src_idx = fdtd_index(g, is + off, js + off, ks + off);
    size_t *pidx = at_xcalloc(n_probes, sizeof(size_t), "probe idx");
    for (size_t p = 0; p < n_probes; p++)
        pidx[p] = fdtd_index(g, ip + off, jp[p] + off, ks + off);

    int rc = 0;
    for (size_t n = 0; n < n_steps; n++) {
        fdtd_update_h_cpml(g, pml);
        fdtd_update_e_cpml(g, pml);
        double t = (double)(n + 1) * g->dt;
        field_of(g, pol)[src_idx] += at_gauss_deriv(t, t0, tau);
        const double *f = field_of(g, pol);
        for (size_t p = 0; p < n_probes; p++) {
            double v = f[pidx[p]];
            if (!isfinite(v)) { rc = 1; }
            hist[p * n_steps + n] = v;
        }
        if (rc) { fprintf(stderr, "non-finite field at step %zu\n", n); break; }
    }
    free(pidx);
    cpml_destroy(pml);
    fdtd_grid_destroy(g);
    return rc;
}

int main(int argc, char **argv) {
    const double dx = 1e-3;
    const int strict = at_has_flag(argc, argv, "--strict");
    const int extended = at_has_flag(argc, argv, "--extended");
    const size_t npml = (size_t)at_get_long(argc, argv, "--npml", 10);
    const long def_stand = extended ? 3 : 6;
    const size_t ds = (size_t)at_get_long(argc, argv, "--ds", def_stand);
    const size_t dp = (size_t)at_get_long(argc, argv, "--dp", def_stand);
    const double gate_db = at_get_double(argc, argv, "--gate-db", -30.0);
    const double max_mb = at_get_double(argc, argv, "--max-mb", 6000.0);
    const char *pol_s = at_get_str(argc, argv, "--pol", "both");
    const long extra_override = at_get_long(argc, argv, "--reference-extra", -1);

    if (npml < 4 || ds < 1 || dp < 1) { fprintf(stderr, "bad --npml/--ds/--dp\n"); return 2; }

    /* angle list */
    double req[MAX_ANGLES]; size_t n_ang = 0;
    if (at_get_str(argc, argv, "--angle", NULL)) {
        req[n_ang++] = at_get_double(argc, argv, "--angle", 0.0);
    } else {
        const double base[] = {0, 15, 30, 45, 60};
        for (size_t i = 0; i < 5; i++) req[n_ang++] = base[i];
        if (extended) { req[n_ang++] = 75; req[n_ang++] = 85; }
    }
    angle_t ang[MAX_ANGLES];
    size_t dy_max = 0;
    for (size_t a = 0; a < n_ang; a++) {
        if (req[a] < 0.0 || req[a] >= 89.0) { fprintf(stderr, "angle must be in [0, 89)\n"); return 2; }
        double dyf = (double)(ds + dp) * tan(req[a] * M_PI / 180.0);
        ang[a].requested_deg = req[a];
        ang[a].dy = (size_t)llround(dyf);
        ang[a].actual_deg = atan((double)ang[a].dy / (double)(ds + dp)) * 180.0 / M_PI;
        if (ang[a].dy > dy_max) dy_max = ang[a].dy;
    }

    /* test-grid geometry */
    const size_t margin = 8;
    const size_t is = npml + margin;               /* source i               */
    const size_t nx = is + ds + npml;              /* +x CPML face at nx-npml */
    const size_t ip = nx - npml - dp;
    const size_t js = npml + margin;
    const size_t ny = js + dy_max + margin + npml;
    const size_t nz = 2 * npml + 2 * margin;
    const size_t ks = nz / 2;
    if (ip <= npml) { fprintf(stderr, "probe inside -x CPML; reduce --dp\n"); return 2; }

    size_t *jp = at_xcalloc(n_ang, sizeof(size_t), "jp");
    for (size_t a = 0; a < n_ang; a++) jp[a] = js + ang[a].dy;

    /* timing */
    fdtd_grid_t *tmp = fdtd_grid_create(4, 4, 4, dx, 0.99);
    if (!tmp) return 1;
    const double dt = tmp->dt;
    fdtd_grid_destroy(tmp);
    const double tau = 8.0 * dt;
    const double t0 = 4.5 * tau;
    const double fc = at_gauss_deriv_peak_freq(tau);
    const double cells_per_step = AT_C0 * dt / dx;
    double path_max = sqrt((double)((ds + dp) * (ds + dp) + dy_max * dy_max));
    size_t n_steps = (size_t)ceil(2.0 * t0 / dt + path_max / cells_per_step) + 60;

    size_t extra = (extra_override >= 0)
        ? (size_t)extra_override
        : (size_t)ceil(0.5 * (double)n_steps * cells_per_step) + 4;

    const size_t rnx = nx + 2 * extra, rny = ny + 2 * extra, rnz = nz + 2 * extra;
    double ref_mb = at_grid_cpml_mb(rnx, rny, rnz, 1);
    double test_mb = at_grid_cpml_mb(nx, ny, nz, 1);

    printf("=== validate_cpml_oblique ===\n");
    printf("dx=%.3g m dt=%.6g s npml=%zu ds=%zu dp=%zu steps=%zu\n", dx, dt, npml, ds, dp, n_steps);
    printf("test grid %zux%zux%zu (%.0f MB), reference %zux%zux%zu (+%zu/side, %.0f MB)\n",
           nx, ny, nz, test_mb, rnx, rny, rnz, extra, ref_mb);
    printf("pulse: differentiated Gaussian, tau=%.4g s, spectral peak %.4g GHz (%.1f cells/wavelength)\n",
           tau, fc / 1e9, AT_C0 / fc / dx);
    if (extra_override >= 0 && (double)(2 * extra) < (double)n_steps * cells_per_step)
        printf("WARNING: reference-extra %zu is smaller than the clean-window requirement; "
               "late-time reference-boundary reflections may enter the window.\n", extra);
    if (ds < 4 || dp < 4)
        printf("NOTE: standoff < 4 cells -- near-field/evanescent coupling into the CPML is included in the result.\n");
    if (ref_mb + test_mb > max_mb) {
        fprintf(stderr, "required memory %.0f MB exceeds --max-mb %.0f\n", ref_mb + test_mb, max_mb);
        return 2;
    }

    int run_te = (strcmp(pol_s, "both") == 0 || strcmp(pol_s, "te") == 0);
    int run_tm = (strcmp(pol_s, "both") == 0 || strcmp(pol_s, "tm") == 0);
    if (!run_te && !run_tm) { fprintf(stderr, "--pol must be te, tm or both\n"); return 2; }

    const char *sum_path = "cpml_oblique_summary.csv";
    FILE *fs = at_xfopen(sum_path);
    fprintf(fs, "pol,requested_deg,actual_deg,dy_cells,npml,ds,dp,max_norm,max_norm_db,rel_l2,spec_db_at_fc,ref_peak,pass\n");

    int n_fail = 0, rc_total = 0;
    double *ht = at_xcalloc(n_ang * n_steps, sizeof(double), "test hist");
    double *hr = at_xcalloc(n_ang * n_steps, sizeof(double), "ref hist");

    for (int pp = 0; pp < 2; pp++) {
        pol_t pol = (pol_t)pp;
        if ((pol == POL_TE && !run_te) || (pol == POL_TM && !run_tm)) continue;
        const char *pname = pol == POL_TE ? "TE" : "TM";
        printf("\n--- %s (%s dipole) ---\n", pname, pol == POL_TE ? "Ez" : "Ey");

        if (run_sim(nx, ny, nz, dx, npml, 0, is, js, ks, ip, jp, n_ang, pol, n_steps, t0, tau, ht) ||
            run_sim(rnx, rny, rnz, dx, npml, extra, is, js, ks, ip, jp, n_ang, pol, n_steps, t0, tau, hr)) {
            rc_total = 1; break;
        }

        char wpath[64];
        snprintf(wpath, sizeof wpath, "cpml_oblique_waveforms_%s.csv", pol == POL_TE ? "te" : "tm");
        FILE *fw = at_xfopen(wpath);
        fprintf(fw, "step,t");
        for (size_t a = 0; a < n_ang; a++)
            fprintf(fw, ",ref_%.0fdeg,test_%.0fdeg", ang[a].requested_deg, ang[a].requested_deg);
        fprintf(fw, "\n");
        for (size_t n = 0; n < n_steps; n++) {
            fprintf(fw, "%zu,%.9e", n, (double)(n + 1) * dt);
            for (size_t a = 0; a < n_ang; a++)
                fprintf(fw, ",%.9e,%.9e", hr[a * n_steps + n], ht[a * n_steps + n]);
            fprintf(fw, "\n");
        }
        if (at_xfclose(fw, wpath)) rc_total = 1;

        printf("%9s %9s %6s %12s %10s %12s %12s  %s\n",
               "req_deg", "act_deg", "dy", "max_norm", "max_dB", "rel_l2", "spec_dB@fc", "result");
        for (size_t a = 0; a < n_ang; a++) {
            const double *r = hr + a * n_steps, *t = ht + a * n_steps;
            double peak = 0, emax = 0, e2 = 0, r2 = 0;
            at_dft_t dr = {0, 0}, de = {0, 0};
            for (size_t n = 0; n < n_steps; n++) {
                double e = t[n] - r[n], tt = (double)(n + 1) * dt;
                if (fabs(r[n]) > peak) peak = fabs(r[n]);
                if (fabs(e) > emax) emax = fabs(e);
                e2 += e * e; r2 += r[n] * r[n];
                at_dft_add(&dr, r[n], tt, fc, dt);
                at_dft_add(&de, e, tt, fc, dt);
            }
            double mn = peak > 0 ? emax / peak : INFINITY;
            double l2 = r2 > 0 ? sqrt(e2 / r2) : INFINITY;
            double sdb = at_db20(sqrt(at_dft_abs2(de) / at_dft_abs2(dr)));
            int pass = isfinite(mn) && at_db20(mn) <= gate_db;
            if (!pass) n_fail++;
            printf("%9.1f %9.2f %6zu %12.4e %10.2f %12.4e %12.2f  %s\n",
                   ang[a].requested_deg, ang[a].actual_deg, ang[a].dy, mn, at_db20(mn), l2, sdb,
                   pass ? "PASS" : "FAIL");
            fprintf(fs, "%s,%.1f,%.4f,%zu,%zu,%zu,%zu,%.9e,%.4f,%.9e,%.4f,%.9e,%d\n",
                    pname, ang[a].requested_deg, ang[a].actual_deg, ang[a].dy, npml, ds, dp,
                    mn, at_db20(mn), l2, sdb, peak, pass);
        }
    }

    free(ht); free(hr); free(jp);
    if (at_xfclose(fs, sum_path)) rc_total = 1;

    if (rc_total) { printf("\nRESULT: ERROR (execution failure)\n"); return 1; }
    printf("\nRESULT: %d probe(s) above the %.1f dB gate%s\n", n_fail, gate_db,
           strict ? "" : " (diagnostic mode: not treated as failure; use --strict)");
    return (strict && n_fail) ? 1 : 0;
}
