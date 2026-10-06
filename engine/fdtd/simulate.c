#include "simulate.h"
#include "grid.h"
#include "cpml.h"
#include "materials.h"
#include "source.h"
#include "probe.h"
#include <stdlib.h>
#include <math.h>

int fdtd_simulate_echo(const fdtd_sim_config_t *cfg, double *out_td, double *sample_rate_out) {
    if (!cfg || !out_td) return -1;
    if (cfg->n_steps == 0) return -1;

    fdtd_grid_t *g = fdtd_grid_create(cfg->nx, cfg->ny, cfg->nz,
                                      cfg->dx, cfg->courant_safety);
    if (!g) return -1;

    if (cfg->src_i > g->nx || cfg->src_j > g->ny || cfg->src_k > g->nz ||
        cfg->probe_i > g->nx || cfg->probe_j > g->ny || cfg->probe_k > g->nz) {
        fdtd_grid_destroy(g);
        return -1;
    }

    cpml_t *pml = cpml_create(g, cpml_default_params(cfg->npml));
    if (!pml) { fdtd_grid_destroy(g); return -1; }

    /* vacuum for now; world.c fills this in phase 2 */
    materials_t *mat = materials_create(g);
    if (!mat) { cpml_destroy(pml); fdtd_grid_destroy(g); return -1; }

    source_gaussian_deriv_params_t sp = { cfg->t0_over_dt * g->dt,
                                          cfg->tau_over_dt * g->dt };
    source_t src = source_make_point(cfg->src_i, cfg->src_j, cfg->src_k, cfg->src_component, SRC_SOFT, 1.0, source_waveform_gaussian_derivative, &sp);

    probe_td_t probe;
    if (probe_td_create(&probe, cfg->probe_i, cfg->probe_j, cfg->probe_k, cfg->probe_component, cfg->n_steps) != 0) {
        materials_destroy(mat);
        cpml_destroy(pml);
        fdtd_grid_destroy(g);
        return -1;
    }

    for (size_t n = 0; n < cfg->n_steps; n++) {
        double t = (double)n * g->dt;
        fdtd_update_h_cpml(g, pml);
        source_apply(&src, g, t);
        fdtd_update_e_cpml_materials(g, pml, mat);
        probe_td_record(&probe, g);
    }

    for (size_t n = 0; n < cfg->n_steps; n++)
        out_td[n] = probe.buffer[n];

    if (sample_rate_out) *sample_rate_out = 1.0 / g->dt;

    probe_td_destroy(&probe);
    materials_destroy(mat);
    cpml_destroy(pml);
    fdtd_grid_destroy(g);
    return 0;
}