#ifndef ENGINE_FDTD_MATERIALS_H
#define ENGINE_FDTD_MATERIALS_H

#include <stddef.h>
#include "grid.h"
#include "cpml.h"

typedef struct materials materials_t;

materials_t *materials_create(const fdtd_grid_t *g);
void materials_destroy(materials_t *mat);

int materials_set_box(materials_t *mat, const fdtd_grid_t *g, size_t i0, size_t i1, size_t j0, size_t j1, size_t k0, size_t k1, double eps_r, double sigma);

void fdtd_update_e_materials(fdtd_grid_t *g, const materials_t *mat);

void fdtd_update_e_cpml_materials(fdtd_grid_t *g, cpml_t *pml, const materials_t *mat);

#endif