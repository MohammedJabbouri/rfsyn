#ifndef ENGINE_FDTD_MATERIALS_H
#define ENGINE_FDTD_MATERIALS_H

#include <stddef.h>
#include <stdint.h>
#include "grid.h"
#include "cpml.h"
#include "field_component.h"

typedef struct materials materials_t;
typedef uint8_t material_id_t;
typedef struct { double Ca, Cb; } material_coeff_t;

materials_t *materials_create(const fdtd_grid_t *g);
void materials_destroy(materials_t *mat);

int materials_register(materials_t *mat, double eps_r, double sigma, material_id_t *id_out);

int materials_assign_box(materials_t *mat, const fdtd_grid_t *g, size_t i0, size_t i1, size_t j0, size_t j1, size_t k0, size_t k1, material_id_t id);
int materials_set_box(materials_t *mat, const fdtd_grid_t *g, size_t i0, size_t i1, size_t j0, size_t j1, size_t k0, size_t k1, double eps_r, double sigma);

int materials_set_component(materials_t *mat, const fdtd_grid_t *g, field_component_t component, size_t i, size_t j, size_t k, double eps_r, double sigma);
int materials_clear_component(materials_t *mat, const fdtd_grid_t *g, field_component_t component, size_t i, size_t j, size_t k);
int materials_get_coefficients(const materials_t *mat, const fdtd_grid_t *g, field_component_t component, size_t i, size_t j, size_t k, material_coeff_t *out);

void fdtd_update_e_materials(fdtd_grid_t *g, const materials_t *mat);
void fdtd_update_e_cpml_materials(fdtd_grid_t *g, cpml_t *pml, const materials_t *mat);

#endif