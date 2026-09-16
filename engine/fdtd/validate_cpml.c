#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "cpml.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static double source_waveform(double t, double t0, double tau) {
    double arg = (t - t0) / tau;
    return -2.0 * arg * exp(-arg * arg);
}

int main(void) {
  const size_t N_REF = 120;
  const size_t N_TEST = 40;
  const double dx = 1e-3;
  const size_t NPML = 8;
  const size_t n_steps = 150;

  size_t src_ref = N_REF / 2, probe_ref_i = N_REF / 2 + 12, probe_ref_j = N_REF/2, probe_ref_k = N_REF/2;
  size_t src_test = N_TEST / 2, probe_test_i = N_TEST / 2 + 12, probe_test_j = N_TEST/2, probe_test_k = N_TEST/2;

  fdtd_grid_t *g_ref = fdtd_grid_create(N_REF, N_REF, N_REF, dx, 0.99);
  fdtd_grid_t *g_test = fdtd_grid_create(N_TEST, N_TEST, N_TEST, dx, 0.99);
  if (!g_ref || !g_test) { fprintf(stderr, "grid alloc failed\n"); return 1; }

  if (fabs(g_ref->dt - g_test->dt) > 1e-30) {
    fprintf(stderr, "dt mismatch between grids | unexpected, check fdtd_grid_create\n");
    return 1;
  }

  cpml_t *pml = cpml_create(g_test, cpml_default_params(NPML));
  if (!pml) { fprintf(stderr, "cpml_create failed\n"); return 1; }

  double t0 = 40.0 * g_ref->dt;
  double tau = 8.0 * g_ref->dt;

  FILE *out = fopen("cpml_validation.csv", "w");
  fprintf(out, "step,t,probe_ref,probe_test,abs_err,db_err\n");

  double max_ref_mag = 0.0;

  for (size_t n = 0; n < n_steps; n++) {
    double t = (double)n * g_ref->dt;
    double src_val = source_waveform(t, t0, tau);

    fdtd_update_h(g_ref);
      fdtd_update_h_cpml(g_test, pml);

    g_ref->ez[fdtd_index(g_ref, src_ref, src_ref, src_ref)] += src_val;
    g_test->ez[fdtd_index(g_test, src_test, src_test, src_test)] += src_val;

    fdtd_update_e(g_ref);
    fdtd_update_e_cpml(g_test, pml);

    double p_ref = g_ref->ez[fdtd_index(g_ref, probe_ref_i, probe_ref_j, probe_ref_k)];
    double p_test = g_test->ez[fdtd_index(g_test, probe_test_i, probe_test_j, probe_test_k)];

    if (fabs(p_ref) > max_ref_mag) max_ref_mag = fabs(p_ref);

    double abs_err = fabs(p_test - p_ref);
    double db_err = (max_ref_mag > 0.0)
    ? 20.0 * log10(abs_err / max_ref_mag + 1e-300) : 0.0;

    fprintf(out, "%zu,%.9e,%.9e,%.9e,%.9e,%.6f\n", n, t, p_ref, p_test, abs_err, db_err);
  }

  fclose(out);

  cpml_destroy(pml);

  fdtd_grid_destroy(g_test);
  fdtd_grid_destroy(g_ref);

  printf("wrote cpml_validation.csv -- inspect db_err column\n");
  printf("expect it to bottom out somewhere near your r0 target (%.0e -> ~%.1f dB)\n", 1e-6, 20.0 * log10(1e-6));

  return 0;
}

