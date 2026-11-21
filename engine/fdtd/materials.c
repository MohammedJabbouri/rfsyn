#include "materials.h"
#include "../core/constants.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MATERIAL_LIMIT 256

typedef struct {
    double eps_r;
    double sigma;
    material_coeff_t coeff;
} material_entry_t;

typedef struct {
    size_t key;
    material_coeff_t coeff;
} override_t;

typedef struct {
    size_t key;
    material_tensor_t coeff;
} tensor_override_t;

struct materials {
    size_t nx;
    size_t ny;
    size_t nz;
    size_t total;
    double dx;
    double dt;
    uint8_t *ids;
    material_entry_t table[MATERIAL_LIMIT];
    size_t table_count;
    override_t *overrides;
    size_t override_count;
    size_t override_capacity;
    tensor_override_t *tensor_overrides;
    size_t tensor_override_count;
    size_t tensor_override_capacity;
};

// Materials must match the grid geometry and time step.
static int grid_matches(const materials_t *mat, const fdtd_grid_t *g) {
    return mat && g &&
           mat->nx == g->nx &&
           mat->ny == g->ny &&
           mat->nz == g->nz &&
           mat->dx == g->dx &&
           mat->dt == g->dt;
}

static int component_number(field_component_t component) {
    switch (component) {
        case FIELD_EX: return 0;
        case FIELD_EY: return 1;
        case FIELD_EZ: return 2;
        default: return -1;
    }
}

// Periodic axes omit the duplicate endpoint from component coordinates.
static int component_valid(const fdtd_grid_t *g, int component,
                           size_t i, size_t j, size_t k) {
    if (!g || component < 0 || component > 2) return 0;

    size_t ni = g->nx;
    size_t nj = g->ny;
    size_t nk = g->nz;

    if (component != 0 && !(g->periodic_axes & 1u)) ni++;
    if (component != 1 && !(g->periodic_axes & 2u)) nj++;
    if (component != 2 && !(g->periodic_axes & 4u)) nk++;

    return i < ni && j < nj && k < nk;
}

static int make_coeff(const materials_t *mat,
                      double eps_r, double sigma,
                      material_coeff_t *out) {
    if (!mat || !out ||
        !isfinite(eps_r) || eps_r < 1.0 ||
        !isfinite(sigma) || sigma < 0.0) {
        return -1;
    }

    double eps = eps_r * EPS0;
    double loss = (sigma * mat->dt) / (2.0 * eps);
    double cb = mat->dt / (eps * mat->dx);

    if (!isfinite(eps) || eps <= 0.0 ||
        !isfinite(loss) ||
        !isfinite(cb) || cb <= 0.0) {
        return -1;
    }

    material_coeff_t coeff = {
        .Ca = (1.0 - loss) / (1.0 + loss),
        .Cb = cb / (1.0 + loss)
    };

    if (!isfinite(coeff.Ca) ||
        !isfinite(coeff.Cb) || coeff.Cb <= 0.0) {
        return -1;
    }

    *out = coeff;
    return 0;
}

// Overrides remain sorted by cell index and electric component.
static size_t lower_bound(const materials_t *mat, size_t key) {
    size_t lo = 0;
    size_t hi = mat->override_count;

    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;

        if (mat->overrides[mid].key < key) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    return lo;
}

static material_coeff_t lookup(const materials_t *mat,
                               size_t idx, int component) {
    if (mat->override_count != 0) {
        size_t key = idx * 3 + (size_t) component;
        size_t pos = lower_bound(mat, key);

        if (pos < mat->override_count &&
            mat->overrides[pos].key == key) {
            return mat->overrides[pos].coeff;
        }
    }

    return mat->table[mat->ids[idx]].coeff;
}

static int reserve_override(materials_t *mat) {
    if (mat->override_count < mat->override_capacity) return 0;

    size_t limit = SIZE_MAX / sizeof(*mat->overrides);
    size_t capacity = mat->override_capacity;

    if (capacity >= limit) return -1;

    size_t next;

    if (capacity == 0) {
        next = limit < 16 ? limit : 16;
    } else {
        next = capacity > limit / 2 ? limit : capacity * 2;
    }

    if (next <= capacity) return -1;

    override_t *grown = realloc(
        mat->overrides,
        next * sizeof(*mat->overrides)
    );

    if (!grown) return -1;

    mat->overrides = grown;
    mat->override_capacity = next;
    return 0;
}

// Tensor overrides share the key/ordering scheme of scalar overrides.
static size_t tensor_lower_bound(const materials_t *mat, size_t key) {
    size_t lo = 0, hi = mat->tensor_override_count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (mat->tensor_overrides[mid].key < key) lo = mid + 1;
        else hi = mid;
    }
    return lo;
}

static const material_tensor_t *tensor_lookup(const materials_t *mat,
                                              size_t idx, int component) {
    if (mat->tensor_override_count == 0) return NULL;
    size_t key = idx * 3 + (size_t) component;
    size_t pos = tensor_lower_bound(mat, key);
    if (pos < mat->tensor_override_count &&
        mat->tensor_overrides[pos].key == key) {
        return &mat->tensor_overrides[pos].coeff;
    }
    return NULL;
}

static int tensor_reserve(materials_t *mat) {
    if (mat->tensor_override_count < mat->tensor_override_capacity) return 0;
    size_t limit = SIZE_MAX / sizeof(*mat->tensor_overrides);
    size_t capacity = mat->tensor_override_capacity;
    if (capacity >= limit) return -1;
    size_t next = capacity == 0
                      ? (limit < 16 ? limit : 16)
                      : (capacity > limit / 2 ? limit : capacity * 2);
    if (next <= capacity) return -1;
    tensor_override_t *grown = realloc(mat->tensor_overrides,
                                       next * sizeof(*mat->tensor_overrides));
    if (!grown) return -1;
    mat->tensor_overrides = grown;
    mat->tensor_override_capacity = next;
    return 0;
}

// Tensor edges need one-cell interior margins on nonperiodic axes so the
// centered cross-curl never samples across a hard boundary.
static int tensor_margin_ok(const fdtd_grid_t *g, int component, size_t i, size_t j, size_t k) {
    (void)component;
    (void)k;
    if (!(g->periodic_axes & 1u) && (i < 1 || i + 2 > g->nx)) return 0;
    if (!(g->periodic_axes & 2u) && (j < 1 || j + 2 > g->ny)) return 0;
    return 1;
}

materials_t *materials_create(const fdtd_grid_t *g) {
    if (!g || g->nx == 0 || g->ny == 0 || g->nz == 0) return NULL;

    if (!isfinite(g->dx) || g->dx <= 0.0 ||
        !isfinite(g->dt) || g->dt <= 0.0) {
        return NULL;
    }

    if (g->nx == SIZE_MAX ||
        g->ny == SIZE_MAX ||
        g->nz == SIZE_MAX) {
        return NULL;
    }

    size_t nx1 = g->nx + 1;
    size_t ny1 = g->ny + 1;
    size_t nz1 = g->nz + 1;

    if (nx1 > SIZE_MAX / ny1) return NULL;

    size_t total = nx1 * ny1;

    if (total > SIZE_MAX / nz1) return NULL;

    total *= nz1;

    // Each padded cell reserves three possible component override keys.
    if (total > SIZE_MAX / 3) return NULL;

    materials_t *mat = calloc(1, sizeof(*mat));
    if (!mat) return NULL;

    mat->nx = g->nx;
    mat->ny = g->ny;
    mat->nz = g->nz;
    mat->total = total;
    mat->dx = g->dx;
    mat->dt = g->dt;

    mat->ids = calloc(total, sizeof(*mat->ids));

    if (!mat->ids) {
        materials_destroy(mat);
        return NULL;
    }

    // Zero-initialized material IDs refer to the first registered material.
    material_id_t vacuum;

    if (materials_register(mat, 1.0, 0.0, &vacuum) != 0) {
        materials_destroy(mat);
        return NULL;
    }

    return mat;
}

void materials_destroy(materials_t *mat) {
    if (!mat) return;
    free(mat->ids);
    free(mat->overrides);
    free(mat->tensor_overrides);
    free(mat);
}

int materials_register(materials_t *mat,
                       double eps_r, double sigma,
                       material_id_t *id_out) {
    material_coeff_t coeff;

    if (!mat || !id_out ||
        make_coeff(mat, eps_r, sigma, &coeff) != 0) {
        return -1;
    }

    for (size_t n = 0; n < mat->table_count; n++) {
        if (mat->table[n].eps_r == eps_r &&
            mat->table[n].sigma == sigma) {
            *id_out = (material_id_t) n;
            return 0;
        }
    }

    if (mat->table_count >= MATERIAL_LIMIT) return -1;

    size_t n = mat->table_count;

    mat->table[n] = (material_entry_t){
        .eps_r = eps_r,
        .sigma = sigma,
        .coeff = coeff
    };

    mat->table_count++;
    *id_out = (material_id_t) n;
    return 0;
}

// AI DEBUG
// Padded box ends are accepted, but periodic starts must be canonical.
static int box_valid(const materials_t *mat, const fdtd_grid_t *g,
                     size_t i0, size_t i1,
                     size_t j0, size_t j1,
                     size_t k0, size_t k1) {
    if (!grid_matches(mat, g)) return 0;

    if (i0 >= i1 || j0 >= j1 || k0 >= k1) return 0;

    if (i1 > g->nx + 1 ||
        j1 > g->ny + 1 ||
        k1 > g->nz + 1) {
        return 0;
    }

    if ((g->periodic_axes & 1u) && i0 >= g->nx) return 0;
    if ((g->periodic_axes & 2u) && j0 >= g->ny) return 0;
    if ((g->periodic_axes & 4u) && k0 >= g->nz) return 0;

    return 1;
}

int materials_assign_box(materials_t *mat, const fdtd_grid_t *g, size_t i0, size_t i1, size_t j0, size_t j1, size_t k0,
                         size_t k1, material_id_t id) {
    if (!box_valid(mat, g, i0, i1, j0, j1, k0, k1) || (size_t) id >= mat->table_count) {
        return -1;
    }

    // Clip periodic endpoints instead of aliasing them onto coordinate zero.
    if ((g->periodic_axes & 1u) && i1 > g->nx) i1 = g->nx;
    if ((g->periodic_axes & 2u) && j1 > g->ny) j1 = g->ny;
    if ((g->periodic_axes & 4u) && k1 > g->nz) k1 = g->nz;

    for (size_t i = i0; i < i1; i++) {
        for (size_t j = j0; j < j1; j++) {
            for (size_t k = k0; k < k1; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                mat->ids[idx] = id;
            }
        }
    }

    // Assignment clears component overrides only inside the clipped box.
    size_t kept = 0;

    for (size_t n = 0; n < mat->override_count; n++) {
        size_t idx = mat->overrides[n].key / 3;
        size_t k = idx % (g->nz + 1);
        size_t row = idx / (g->nz + 1);
        size_t j = row % (g->ny + 1);
        size_t i = row / (g->ny + 1);

        if (i >= i0 && i < i1 && j >= j0 && j < j1 && k >= k0 && k < k1) {
            continue;
        }

        mat->overrides[kept++] = mat->overrides[n];
    }

    mat->override_count = kept;

    // Assignment clears tensor overrides only inside the clipped box.
    size_t tensor_kept = 0;

    for (size_t n = 0; n < mat->tensor_override_count; n++) {
        size_t idx = mat->tensor_overrides[n].key / 3;
        size_t k = idx % (g->nz + 1);
        size_t row = idx / (g->nz + 1);
        size_t j = row % (g->ny + 1);
        size_t i = row / (g->ny + 1);

        if (i >= i0 && i < i1 && j >= j0 && j < j1 && k >= k0 && k < k1) {
            continue;
        }

        mat->tensor_overrides[tensor_kept++] = mat->tensor_overrides[n];
    }

    mat->tensor_override_count = tensor_kept;
    return 0;
}

int materials_set_box(materials_t *mat, const fdtd_grid_t *g,
                      size_t i0, size_t i1,
                      size_t j0, size_t j1,
                      size_t k0, size_t k1,
                      double eps_r, double sigma) {
    if (!box_valid(mat, g, i0, i1, j0, j1, k0, k1)) return -1;

    material_id_t id;

    if (materials_register(mat, eps_r, sigma, &id) != 0) return -1;

    return materials_assign_box(
        mat, g,
        i0, i1,
        j0, j1,
        k0, k1,
        id
    );
}

int materials_set_component(materials_t *mat, const fdtd_grid_t *g,
                            field_component_t component,
                            size_t i, size_t j, size_t k,
                            double eps_r, double sigma) {
    int c = component_number(component);
    material_coeff_t coeff;

    if (!grid_matches(mat, g) ||
        !component_valid(g, c, i, j, k) ||
        make_coeff(mat, eps_r, sigma, &coeff) != 0) {
        return -1;
    }

    size_t key = fdtd_index(g, i, j, k) * 3 + (size_t) c;
    size_t pos = lower_bound(mat, key);

    if (pos < mat->override_count &&
        mat->overrides[pos].key == key) {
        mat->overrides[pos].coeff = coeff;
        return 0;
    }

    if (reserve_override(mat) != 0) return -1;

    size_t trailing = mat->override_count - pos;

    if (trailing != 0) {
        memmove(
            mat->overrides + pos + 1,
            mat->overrides + pos,
            trailing * sizeof(*mat->overrides)
        );
    }

    mat->overrides[pos] = (override_t){
        .key = key,
        .coeff = coeff
    };

    mat->override_count++;
    return 0;
}

int materials_clear_component(materials_t *mat, const fdtd_grid_t *g,
                              field_component_t component,
                              size_t i, size_t j, size_t k) {
    int c = component_number(component);

    if (!grid_matches(mat, g) ||
        !component_valid(g, c, i, j, k)) {
        return -1;
    }

    size_t key = fdtd_index(g, i, j, k) * 3 + (size_t) c;
    size_t pos = lower_bound(mat, key);

    if (pos < mat->override_count &&
        mat->overrides[pos].key == key) {
        size_t trailing = mat->override_count - pos - 1;

        if (trailing != 0) {
            memmove(
                mat->overrides + pos,
                mat->overrides + pos + 1,
                trailing * sizeof(*mat->overrides)
            );
        }

        mat->override_count--;
    }

    return 0;
}

int materials_set_component_tensor(materials_t *mat, const fdtd_grid_t *g,
                                   field_component_t component, size_t i, size_t j, size_t k,
                                   double m11, double m12, double m22) {
    int c = component_number(component);
    if (c == 2) return -1; /* Ez tensor smoothing is out of scope (v1) */

    if (!grid_matches(mat, g) ||
        !component_valid(g, c, i, j, k) ||
        !tensor_margin_ok(g, c, i, j, k) ||
        !isfinite(m11) || !isfinite(m12) || !isfinite(m22) ||
        !(m11 > 0.0) || !(m22 > 0.0) ||
        !(m11 * m22 - m12 * m12 > 0.0)) {
        /* SPD required */
        return -1;
    }

    double scale = mat->dt / (EPS0 * mat->dx);
    if (!isfinite(scale) || scale <= 0.0) return -1;

    size_t key = fdtd_index(g, i, j, k) * 3 + (size_t) c;
    size_t pos = tensor_lower_bound(mat, key);

    material_tensor_t coeff = {
        .m11 = m11, .m12 = m12, .m22 = m22,
        .scale = scale
    };

    if (pos < mat->tensor_override_count &&
        mat->tensor_overrides[pos].key == key) {
        mat->tensor_overrides[pos].coeff = coeff;
        return 0;
    }
    if (tensor_reserve(mat) != 0) return -1;

    size_t trailing = mat->tensor_override_count - pos;
    if (trailing != 0) {
        memmove(mat->tensor_overrides + pos + 1,
                mat->tensor_overrides + pos,
                trailing * sizeof(*mat->tensor_overrides));
    }
    mat->tensor_overrides[pos] = (tensor_override_t)
    {
        .key = key, .coeff = coeff
    };
    mat->tensor_override_count++;
    return 0;
}

int materials_clear_component_tensor(materials_t *mat, const fdtd_grid_t *g,
                                     field_component_t component, size_t i, size_t j, size_t k) {
    int c = component_number(component);
    if (!grid_matches(mat, g) || !component_valid(g, c, i, j, k)) return -1;
    size_t key = fdtd_index(g, i, j, k) * 3 + (size_t) c;
    size_t pos = tensor_lower_bound(mat, key);
    if (pos < mat->tensor_override_count &&
        mat->tensor_overrides[pos].key == key) {
        size_t trailing = mat->tensor_override_count - pos - 1;
        if (trailing != 0) {
            memmove(mat->tensor_overrides + pos,
                    mat->tensor_overrides + pos + 1,
                    trailing * sizeof(*mat->tensor_overrides));
        }
        mat->tensor_override_count--;
    }
    return 0;
}

int materials_get_tensor(const materials_t *mat, const fdtd_grid_t *g,
                         field_component_t component, size_t i, size_t j, size_t k,
                         material_tensor_t *out) {
    int c = component_number(component);
    if (!out || !grid_matches(mat, g) || !component_valid(g, c, i, j, k))
        return -1;
    const material_tensor_t *t = tensor_lookup(mat, fdtd_index(g, i, j, k), c);
    if (!t) return 1; /* no tensor at this edge */
    *out = *t;
    return 0;
}

int materials_get_coefficients(const materials_t *mat,
                               const fdtd_grid_t *g,
                               field_component_t component,
                               size_t i, size_t j, size_t k,
                               material_coeff_t *out) {
    int c = component_number(component);

    if (!out ||
        !grid_matches(mat, g) ||
        !component_valid(g, c, i, j, k)) {
        return -1;
    }

    *out = lookup(mat, fdtd_index(g, i, j, k), c);
    return 0;
}

// Update loops exclude zero on nonperiodic backward-difference axes.
static inline size_t material_previous(size_t q, size_t extent) {
    return q == 0 ? extent - 1 : q - 1;
}

#define EX(i, j, k) g->ex[fdtd_index(g, (i), (j), (k))]
#define EY(i, j, k) g->ey[fdtd_index(g, (i), (j), (k))]
#define EZ(i, j, k) g->ez[fdtd_index(g, (i), (j), (k))]
#define HX(i, j, k) g->hx[fdtd_index(g, (i), (j), (k))]
#define HY(i, j, k) g->hy[fdtd_index(g, (i), (j), (k))]
#define HZ(i, j, k) g->hz[fdtd_index(g, (i), (j), (k))]

void fdtd_update_e_materials(fdtd_grid_t *g, const materials_t *mat) {
    if (!grid_matches(mat, g)) return;

    size_t i0 = (g->periodic_axes & 1u) ? 0u : 1u;
    size_t j0 = (g->periodic_axes & 2u) ? 0u : 1u;
    size_t k0 = (g->periodic_axes & 4u) ? 0u : 1u;

    for (size_t i = 0; i < g->nx; i++) {
        for (size_t j = j0; j < g->ny; j++) {
            size_t jm = material_previous(j, g->ny);

            for (size_t k = k0; k < g->nz; k++) {
                size_t km = material_previous(k, g->nz);
                size_t idx = fdtd_index(g, i, j, k);
                material_coeff_t coeff = lookup(mat, idx, 0);

                double termA = HZ(i, j, k) - HZ(i, jm, k);
                double termB = HY(i, j, k) - HY(i, j, km);

                const material_tensor_t *t = tensor_lookup(mat, idx, 0);
                if (t) {
                    size_t ip = (i + 1 == g->nx) ? 0 : i + 1; /* periodic wrap */
                    size_t im2 = material_previous(i, g->nx);
                    /* centered dx*dHz/dx at the Ex node */
                    double dcross = 0.25 * (HZ(ip, j, k) + HZ(ip, jm, k) - HZ(im2, j, k) - HZ(im2, jm, k));
                    EX(i, j, k) = EX(i, j, k) + t->scale * (t->m11 * termA - t->m12 * dcross);
                } else {
                    EX(i, j, k) = coeff.Ca * EX(i, j, k) + coeff.Cb * (termA - termB);
                }
            }
        }
    }

    for (size_t i = i0; i < g->nx; i++) {
        size_t im = material_previous(i, g->nx);

        for (size_t j = 0; j < g->ny; j++) {
            for (size_t k = k0; k < g->nz; k++) {
                size_t km = material_previous(k, g->nz);
                size_t idx = fdtd_index(g, i, j, k);
                material_coeff_t coeff = lookup(mat, idx, 1);

                double termA = HX(i, j, k) - HX(i, j, km);
                double termB = HZ(i, j, k) - HZ(im, j, k);

                const material_tensor_t *t = tensor_lookup(mat, idx, 1);
                if (t) {
                    size_t jp = (j + 1 == g->ny) ? 0 : j + 1;
                    size_t jm2 = material_previous(j, g->ny);
                    size_t jmm = material_previous(jm2, g->ny);
                    /* centered dx*dHz/dy at the Ey node */
                    double dcross = 0.25 * (HZ(i, jp, k) + HZ(im, jp, k) - HZ(i, jmm, k) - HZ(im, jmm, k));
                    EY(i, j, k) = EY(i, j, k) + t->scale * (t->m11 * (-termB) + t->m12 * dcross);
                } else {
                    EY(i, j, k) = coeff.Ca * EY(i, j, k) + coeff.Cb * (termA - termB);
                }
            }
        }
    }

    for (size_t i = i0; i < g->nx; i++) {
        size_t im = material_previous(i, g->nx);

        for (size_t j = j0; j < g->ny; j++) {
            size_t jm = material_previous(j, g->ny);

            for (size_t k = 0; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                material_coeff_t coeff = lookup(mat, idx, 2);

                double termA = HY(i, j, k) - HY(im, j, k);
                double termB = HX(i, j, k) - HX(i, jm, k);

                EZ(i, j, k) = coeff.Ca * EZ(i, j, k) + coeff.Cb * (termA - termB);
            }
        }
    }
}

void fdtd_update_e_cpml_materials(fdtd_grid_t *g, cpml_t *pml,
                                  const materials_t *mat) {
    if (!grid_matches(mat, g) || !pml || pml->dx != g->dx) return;

    size_t i0 = (g->periodic_axes & 1u) ? 0u : 1u;
    size_t j0 = (g->periodic_axes & 2u) ? 0u : 1u;
    size_t k0 = (g->periodic_axes & 4u) ? 0u : 1u;

    for (size_t i = 0; i < g->nx; i++) {
        for (size_t j = j0; j < g->ny; j++) {
            size_t jm = material_previous(j, g->ny);

            for (size_t k = k0; k < g->nz; k++) {
                size_t km = material_previous(k, g->nz);
                size_t idx = fdtd_index(g, i, j, k);
                material_coeff_t coeff = lookup(mat, idx, 0);

                double termA = HZ(i, j, k) - HZ(i, jm, k);
                double termB = HY(i, j, k) - HY(i, j, km);

                termA = cpml_term(
                    termA,
                    pml->kappa_inv_m_y[j],
                    pml->b_m_y[j],
                    pml->a_m_y[j],
                    g->dx,
                    &pml->psi_ex_y[idx]
                );

                termB = cpml_term(
                    termB,
                    pml->kappa_inv_m_z[k],
                    pml->b_m_z[k],
                    pml->a_m_z[k],
                    g->dx,
                    &pml->psi_ex_z[idx]
                );

                EX(i, j, k) = coeff.Ca * EX(i, j, k) +
                              coeff.Cb * (termA - termB);
            }
        }
    }

    for (size_t i = i0; i < g->nx; i++) {
        size_t im = material_previous(i, g->nx);

        for (size_t j = 0; j < g->ny; j++) {
            for (size_t k = k0; k < g->nz; k++) {
                size_t km = material_previous(k, g->nz);
                size_t idx = fdtd_index(g, i, j, k);
                material_coeff_t coeff = lookup(mat, idx, 1);

                double termA = HX(i, j, k) - HX(i, j, km);
                double termB = HZ(i, j, k) - HZ(im, j, k);

                termA = cpml_term(
                    termA,
                    pml->kappa_inv_m_z[k],
                    pml->b_m_z[k],
                    pml->a_m_z[k],
                    g->dx,
                    &pml->psi_ey_z[idx]
                );

                termB = cpml_term(
                    termB,
                    pml->kappa_inv_m_x[i],
                    pml->b_m_x[i],
                    pml->a_m_x[i],
                    g->dx,
                    &pml->psi_ey_x[idx]
                );

                EY(i, j, k) = coeff.Ca * EY(i, j, k) +
                              coeff.Cb * (termA - termB);
            }
        }
    }

    for (size_t i = i0; i < g->nx; i++) {
        size_t im = material_previous(i, g->nx);

        for (size_t j = j0; j < g->ny; j++) {
            size_t jm = material_previous(j, g->ny);

            for (size_t k = 0; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                material_coeff_t coeff = lookup(mat, idx, 2);

                double termA = HY(i, j, k) - HY(im, j, k);
                double termB = HX(i, j, k) - HX(i, jm, k);

                termA = cpml_term(
                    termA,
                    pml->kappa_inv_m_x[i],
                    pml->b_m_x[i],
                    pml->a_m_x[i],
                    g->dx,
                    &pml->psi_ez_x[idx]
                );

                termB = cpml_term(
                    termB,
                    pml->kappa_inv_m_y[j],
                    pml->b_m_y[j],
                    pml->a_m_y[j],
                    g->dx,
                    &pml->psi_ez_y[idx]
                );

                EZ(i, j, k) = coeff.Ca * EZ(i, j, k) +
                              coeff.Cb * (termA - termB);
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
