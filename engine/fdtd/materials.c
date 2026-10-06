#include "materials.h"
#include "../core/constants.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MATERIAL_LIMIT 256

typedef struct { double eps_r, sigma; material_coeff_t coeff; } material_entry_t;
typedef struct { size_t key; material_coeff_t coeff; } override_t;

struct materials {
    size_t nx, ny, nz, total;
    double dx, dt;
    uint8_t *ids;
    material_entry_t table[MATERIAL_LIMIT];
    size_t table_count;
    override_t *overrides;
    size_t override_count, override_capacity;
};

static int grid_matches(const materials_t *mat, const fdtd_grid_t *g) {
    return mat && g && mat->nx == g->nx && mat->ny == g->ny && mat->nz == g->nz && mat->dx == g->dx && mat->dt == g->dt;
}

static int component_number(field_component_t component) {
    switch (component) {
        case FIELD_EX: return 0;
        case FIELD_EY: return 1;
        case FIELD_EZ: return 2;
        default: return -1;
    }
}

static int component_valid(const fdtd_grid_t *g, int c, size_t i, size_t j, size_t k) {
    if (!g || c < 0 || c > 2) return 0;
    if (c == 0) return i < g->nx && j <= g->ny && k <= g->nz;
    if (c == 1) return i <= g->nx && j < g->ny && k <= g->nz;
    return i <= g->nx && j <= g->ny && k < g->nz;
}

static int make_coeff(const materials_t *mat, double eps_r, double sigma, material_coeff_t *out) {
    if (!mat || !out || !isfinite(eps_r) || eps_r < 1.0 || !isfinite(sigma) || sigma < 0.0) return -1;
    double eps = eps_r * EPS0;
    double loss = (sigma * mat->dt) / (2.0 * eps);
    double cb = mat->dt / (eps * mat->dx);
    if (!isfinite(eps) || eps <= 0.0 || !isfinite(loss) || !isfinite(cb) || cb <= 0.0) return -1;
    out->Ca = (1.0 - loss) / (1.0 + loss);
    out->Cb = cb / (1.0 + loss);
    return isfinite(out->Ca) && isfinite(out->Cb) && out->Cb > 0.0 ? 0 : -1;
}

static size_t lower_bound(const materials_t *mat, size_t key) {
    size_t lo = 0, hi = mat->override_count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (mat->overrides[mid].key < key) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

static material_coeff_t lookup(const materials_t *mat, size_t idx, int component) {
    if (mat->override_count) {
        size_t key = idx * 3 + (size_t)component;
        size_t pos = lower_bound(mat, key);
        if (pos < mat->override_count && mat->overrides[pos].key == key) return mat->overrides[pos].coeff;
    }
    return mat->table[mat->ids[idx]].coeff;
}

materials_t *materials_create(const fdtd_grid_t *g) {
    if (!g || g->nx == 0 || g->ny == 0 || g->nz == 0) return NULL;
    if (!isfinite(g->dx) || g->dx <= 0.0 || !isfinite(g->dt) || g->dt <= 0.0) return NULL;
    if (g->nx == SIZE_MAX || g->ny == SIZE_MAX || g->nz == SIZE_MAX) return NULL;
    size_t nx1 = g->nx + 1, ny1 = g->ny + 1, nz1 = g->nz + 1;
    if (nx1 > SIZE_MAX / ny1) return NULL;
    size_t total = nx1 * ny1;
    if (total > SIZE_MAX / nz1) return NULL;
    total *= nz1;
    if (total > SIZE_MAX / 3) return NULL;

    materials_t *mat = calloc(1, sizeof(*mat));
    if (!mat) return NULL;
    mat->nx = g->nx; mat->ny = g->ny; mat->nz = g->nz;
    mat->dx = g->dx; mat->dt = g->dt; mat->total = total;
    mat->ids = calloc(total, sizeof(*mat->ids));
    if (!mat->ids) { free(mat); return NULL; }
    material_id_t vacuum;
    if (materials_register(mat, 1.0, 0.0, &vacuum) != 0) { materials_destroy(mat); return NULL; }
    return mat;
}

void materials_destroy(materials_t *mat) {
    if (!mat) return;
    free(mat->ids);
    free(mat->overrides);
    free(mat);
}

int materials_register(materials_t *mat, double eps_r, double sigma, material_id_t *id_out) {
    material_coeff_t coeff;
    if (!mat || !id_out || make_coeff(mat, eps_r, sigma, &coeff) != 0) return -1;
    for (size_t n = 0; n < mat->table_count; n++) {
        if (mat->table[n].eps_r == eps_r && mat->table[n].sigma == sigma) { *id_out = (material_id_t)n; return 0; }
    }
    if (mat->table_count == MATERIAL_LIMIT) return -1;
    size_t n = mat->table_count++;
    mat->table[n] = (material_entry_t){eps_r, sigma, coeff};
    *id_out = (material_id_t)n;
    return 0;
}

static int box_valid(const materials_t *mat, const fdtd_grid_t *g, size_t i0, size_t i1, size_t j0, size_t j1, size_t k0, size_t k1) {
    return grid_matches(mat, g) && i0 < i1 && j0 < j1 && k0 < k1 && i1 <= g->nx + 1 && j1 <= g->ny + 1 && k1 <= g->nz + 1;
}

int materials_assign_box(materials_t *mat, const fdtd_grid_t *g, size_t i0, size_t i1, size_t j0, size_t j1, size_t k0, size_t k1, material_id_t id) {
    if (!box_valid(mat, g, i0, i1, j0, j1, k0, k1) || (size_t)id >= mat->table_count) return -1;
    for (size_t i = i0; i < i1; i++) {
        for (size_t j = j0; j < j1; j++) {
            size_t idx = fdtd_index(g, i, j, k0);
            memset(mat->ids + idx, id, k1 - k0);
        }
    }
    size_t kept = 0;
    for (size_t n = 0; n < mat->override_count; n++) {
        size_t idx = mat->overrides[n].key / 3;
        size_t k = idx % (g->nz + 1);
        size_t row = idx / (g->nz + 1);
        size_t j = row % (g->ny + 1);
        size_t i = row / (g->ny + 1);
        if (i >= i0 && i < i1 && j >= j0 && j < j1 && k >= k0 && k < k1) continue;
        mat->overrides[kept++] = mat->overrides[n];
    }
    mat->override_count = kept;
    return 0;
}

int materials_set_box(materials_t *mat, const fdtd_grid_t *g, size_t i0, size_t i1, size_t j0, size_t j1, size_t k0, size_t k1, double eps_r, double sigma) {
    if (!box_valid(mat, g, i0, i1, j0, j1, k0, k1)) return -1;
    material_id_t id;
    if (materials_register(mat, eps_r, sigma, &id) != 0) return -1;
    return materials_assign_box(mat, g, i0, i1, j0, j1, k0, k1, id);
}

int materials_set_component(materials_t *mat, const fdtd_grid_t *g, field_component_t component, size_t i, size_t j, size_t k, double eps_r, double sigma) {
    int c = component_number(component);
    material_coeff_t coeff;
    if (!grid_matches(mat, g) || !component_valid(g, c, i, j, k) || make_coeff(mat, eps_r, sigma, &coeff) != 0) return -1;
    size_t key = fdtd_index(g, i, j, k) * 3 + (size_t)c;
    size_t pos = lower_bound(mat, key);
    if (pos < mat->override_count && mat->overrides[pos].key == key) { mat->overrides[pos].coeff = coeff; return 0; }
    if (mat->override_count == mat->override_capacity) {
        size_t limit = SIZE_MAX / sizeof(*mat->overrides);
        size_t cap = mat->override_capacity;
        if (cap == limit) return -1;
        size_t next = cap ? (cap > limit / 2 ? limit : cap * 2) : 16;
        if (next > limit) return -1;
        override_t *grown = realloc(mat->overrides, next * sizeof(*grown));
        if (!grown) return -1;
        mat->overrides = grown;
        mat->override_capacity = next;
    }
    memmove(mat->overrides + pos + 1, mat->overrides + pos, (mat->override_count - pos) * sizeof(*mat->overrides));
    mat->overrides[pos] = (override_t){key, coeff};
    mat->override_count++;
    return 0;
}

int materials_clear_component(materials_t *mat, const fdtd_grid_t *g, field_component_t component, size_t i, size_t j, size_t k) {
    int c = component_number(component);
    if (!grid_matches(mat, g) || !component_valid(g, c, i, j, k)) return -1;
    size_t key = fdtd_index(g, i, j, k) * 3 + (size_t)c;
    size_t pos = lower_bound(mat, key);
    if (pos < mat->override_count && mat->overrides[pos].key == key) {
        memmove(mat->overrides + pos, mat->overrides + pos + 1, (mat->override_count - pos - 1) * sizeof(*mat->overrides));
        mat->override_count--;
    }
    return 0;
}

int materials_get_coefficients(const materials_t *mat, const fdtd_grid_t *g, field_component_t component, size_t i, size_t j, size_t k, material_coeff_t *out) {
    int c = component_number(component);
    if (!out || !grid_matches(mat, g) || !component_valid(g, c, i, j, k)) return -1;
    *out = lookup(mat, fdtd_index(g, i, j, k), c);
    return 0;
}

#define EX(i,j,k) g->ex[fdtd_index(g,(i),(j),(k))]
#define EY(i,j,k) g->ey[fdtd_index(g,(i),(j),(k))]
#define EZ(i,j,k) g->ez[fdtd_index(g,(i),(j),(k))]
#define HX(i,j,k) g->hx[fdtd_index(g,(i),(j),(k))]
#define HY(i,j,k) g->hy[fdtd_index(g,(i),(j),(k))]
#define HZ(i,j,k) g->hz[fdtd_index(g,(i),(j),(k))]

void fdtd_update_e_materials(fdtd_grid_t *g, const materials_t *mat) {
    if (!grid_matches(mat, g)) return;
    for (size_t i = 0; i < g->nx; i++) {
        for (size_t j = 1; j < g->ny; j++) {
            for (size_t k = 1; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                material_coeff_t coeff = lookup(mat, idx, 0);
                double termA = HZ(i,j,k) - HZ(i,j-1,k);
                double termB = HY(i,j,k) - HY(i,j,k-1);
                EX(i,j,k) = coeff.Ca * EX(i,j,k) + coeff.Cb * (termA - termB);
            }
        }
    }
    for (size_t i = 1; i < g->nx; i++) {
        for (size_t j = 0; j < g->ny; j++) {
            for (size_t k = 1; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                material_coeff_t coeff = lookup(mat, idx, 1);
                double termA = HX(i,j,k) - HX(i,j,k-1);
                double termB = HZ(i,j,k) - HZ(i-1,j,k);
                EY(i,j,k) = coeff.Ca * EY(i,j,k) + coeff.Cb * (termA - termB);
            }
        }
    }
    for (size_t i = 1; i < g->nx; i++) {
        for (size_t j = 1; j < g->ny; j++) {
            for (size_t k = 0; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                material_coeff_t coeff = lookup(mat, idx, 2);
                double termA = HY(i,j,k) - HY(i-1,j,k);
                double termB = HX(i,j,k) - HX(i,j-1,k);
                EZ(i,j,k) = coeff.Ca * EZ(i,j,k) + coeff.Cb * (termA - termB);
            }
        }
    }
}

void fdtd_update_e_cpml_materials(fdtd_grid_t *g, cpml_t *pml, const materials_t *mat) {
    if (!grid_matches(mat, g) || !pml || pml->dx != g->dx) return;
    for (size_t i = 0; i < g->nx; i++) {
        for (size_t j = 1; j < g->ny; j++) {
            for (size_t k = 1; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                material_coeff_t coeff = lookup(mat, idx, 0);
                double termA = HZ(i,j,k) - HZ(i,j-1,k);
                double termB = HY(i,j,k) - HY(i,j,k-1);
                termA = cpml_term(termA, pml->kappa_inv_m_y[j], pml->b_m_y[j], pml->a_m_y[j], g->dx, &pml->psi_ex_y[idx]);
                termB = cpml_term(termB, pml->kappa_inv_m_z[k], pml->b_m_z[k], pml->a_m_z[k], g->dx, &pml->psi_ex_z[idx]);
                EX(i,j,k) = coeff.Ca * EX(i,j,k) + coeff.Cb * (termA - termB);
            }
        }
    }
    for (size_t i = 1; i < g->nx; i++) {
        for (size_t j = 0; j < g->ny; j++) {
            for (size_t k = 1; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                material_coeff_t coeff = lookup(mat, idx, 1);
                double termA = HX(i,j,k) - HX(i,j,k-1);
                double termB = HZ(i,j,k) - HZ(i-1,j,k);
                termA = cpml_term(termA, pml->kappa_inv_m_z[k], pml->b_m_z[k], pml->a_m_z[k], g->dx, &pml->psi_ey_z[idx]);
                termB = cpml_term(termB, pml->kappa_inv_m_x[i], pml->b_m_x[i], pml->a_m_x[i], g->dx, &pml->psi_ey_x[idx]);
                EY(i,j,k) = coeff.Ca * EY(i,j,k) + coeff.Cb * (termA - termB);
            }
        }
    }
    for (size_t i = 1; i < g->nx; i++) {
        for (size_t j = 1; j < g->ny; j++) {
            for (size_t k = 0; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                material_coeff_t coeff = lookup(mat, idx, 2);
                double termA = HY(i,j,k) - HY(i-1,j,k);
                double termB = HX(i,j,k) - HX(i,j-1,k);
                termA = cpml_term(termA, pml->kappa_inv_m_x[i], pml->b_m_x[i], pml->a_m_x[i], g->dx, &pml->psi_ez_x[idx]);
                termB = cpml_term(termB, pml->kappa_inv_m_y[j], pml->b_m_y[j], pml->a_m_y[j], g->dx, &pml->psi_ez_y[idx]);
                EZ(i,j,k) = coeff.Ca * EZ(i,j,k) + coeff.Cb * (termA - termB);
            }
        }
    }
}

#undef EX
#undef EY
#undef EZ
#undef HX
#undef HY
#undef HZ