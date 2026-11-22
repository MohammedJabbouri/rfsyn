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

/* Materials must match the grid geometry and time step. */
static int grid_matches(
    const materials_t *mat,
    const fdtd_grid_t *g) {
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

/* Periodic axes omit duplicate endpoint component coordinates. */
static int component_valid(
    const fdtd_grid_t *g,
    int component,
    size_t i,
    size_t j,
    size_t k) {
    if (!g || component < 0 || component > 2) {
        return 0;
    }

    size_t ni = g->nx;
    size_t nj = g->ny;
    size_t nk = g->nz;

    if (component != 0 && !(g->periodic_axes & 1u)) ni++;
    if (component != 1 && !(g->periodic_axes & 2u)) nj++;
    if (component != 2 && !(g->periodic_axes & 4u)) nk++;

    return i < ni && j < nj && k < nk;
}

static int make_coeff(
    const materials_t *mat,
    double eps_r,
    double sigma,
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

/* Scalar overrides are sorted by padded cell index, then component. */
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

static size_t tensor_lower_bound(
    const materials_t *mat,
    size_t key) {
    size_t lo = 0;
    size_t hi = mat->tensor_override_count;

    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;

        if (mat->tensor_overrides[mid].key < key) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    return lo;
}

/* Returns the scalar coefficient, independently of tensor overrides. */
static material_coeff_t lookup(
    const materials_t *mat,
    size_t idx,
    int component) {
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

static const material_tensor_t *tensor_lookup(
    const materials_t *mat,
    size_t idx,
    int component) {
    if (mat->tensor_override_count == 0) {
        return NULL;
    }

    size_t key = idx * 3 + (size_t) component;
    size_t pos = tensor_lower_bound(mat, key);

    if (pos < mat->tensor_override_count &&
        mat->tensor_overrides[pos].key == key) {
        return &mat->tensor_overrides[pos].coeff;
    }

    return NULL;
}

static int reserve_override(materials_t *mat) {
    if (mat->override_count < mat->override_capacity) {
        return 0;
    }

    size_t limit = SIZE_MAX / sizeof(*mat->overrides);
    size_t capacity = mat->override_capacity;
    size_t next;

    if (capacity >= limit) {
        return -1;
    }

    if (capacity == 0) {
        next = limit < 16 ? limit : 16;
    } else {
        next = capacity > limit / 2 ? limit : capacity * 2;
    }

    if (next <= capacity) {
        return -1;
    }

    override_t *grown = realloc(
        mat->overrides,
        next * sizeof(*mat->overrides));

    if (!grown) {
        return -1;
    }

    mat->overrides = grown;
    mat->override_capacity = next;
    return 0;
}

static int tensor_reserve(materials_t *mat) {
    if (mat->tensor_override_count < mat->tensor_override_capacity) {
        return 0;
    }

    size_t limit = SIZE_MAX / sizeof(*mat->tensor_overrides);
    size_t capacity = mat->tensor_override_capacity;
    size_t next;

    if (capacity >= limit) {
        return -1;
    }

    if (capacity == 0) {
        next = limit < 16 ? limit : 16;
    } else {
        next = capacity > limit / 2 ? limit : capacity * 2;
    }

    if (next <= capacity) {
        return -1;
    }

    tensor_override_t *grown = realloc(
        mat->tensor_overrides,
        next * sizeof(*mat->tensor_overrides));

    if (!grown) {
        return -1;
    }

    mat->tensor_overrides = grown;
    mat->tensor_override_capacity = next;
    return 0;
}

/*
 * Conservative x/y interior margins for the tensor interpolation.
 * No z-neighbor samples occur in the restricted tensor stencil.
 */
static int tensor_margin_ok(
    const fdtd_grid_t *g,
    size_t i,
    size_t j) {
    if (!(g->periodic_axes & 1u)) {
        if (g->nx < 3 || i == 0 || i >= g->nx - 1) {
            return 0;
        }
    }

    if (!(g->periodic_axes & 2u)) {
        if (g->ny < 3 || j == 0 || j >= g->ny - 1) {
            return 0;
        }
    }

    return 1;
}

/*
 * Test positive definiteness without forming a*c - b*b.
 * Extremely ill-conditioned inputs may be conservatively rejected
 * because this is a floating-point test.
 */
static int tensor_spd(double a, double b, double c) {
    if (!isfinite(a) || !isfinite(b) || !isfinite(c) ||
        a <= 0.0 || c <= 0.0) {
        return 0;
    }

    return fabs(b) / sqrt(a) < sqrt(c);
}

materials_t *materials_create(const fdtd_grid_t *g) {
    if (!g || g->nx == 0 || g->ny == 0 || g->nz == 0) {
        return NULL;
    }

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

    if (nx1 > SIZE_MAX / ny1) {
        return NULL;
    }

    size_t total = nx1 * ny1;

    if (total > SIZE_MAX / nz1) {
        return NULL;
    }

    total *= nz1;

    /* Reserve three possible component keys per padded cell. */
    if (total > SIZE_MAX / 3) {
        return NULL;
    }

    materials_t *mat = calloc(1, sizeof(*mat));

    if (!mat) {
        return NULL;
    }

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

    material_id_t vacuum;

    if (materials_register(mat, 1.0, 0.0, &vacuum) != 0) {
        materials_destroy(mat);
        return NULL;
    }

    return mat;
}

void materials_destroy(materials_t *mat) {
    if (!mat) {
        return;
    }

    free(mat->ids);
    free(mat->overrides);
    free(mat->tensor_overrides);
    free(mat);
}

int materials_register(
    materials_t *mat,
    double eps_r,
    double sigma,
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

    if (mat->table_count >= MATERIAL_LIMIT) {
        return -1;
    }

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

/* Padded box ends are accepted; periodic starts must be canonical. */
static int box_valid(
    const materials_t *mat,
    const fdtd_grid_t *g,
    size_t i0,
    size_t i1,
    size_t j0,
    size_t j1,
    size_t k0,
    size_t k1) {
    if (!grid_matches(mat, g)) {
        return 0;
    }

    if (i0 >= i1 || j0 >= j1 || k0 >= k1) {
        return 0;
    }

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

static int key_inside_box(
    const fdtd_grid_t *g,
    size_t key,
    size_t i0,
    size_t i1,
    size_t j0,
    size_t j1,
    size_t k0,
    size_t k1) {
    size_t idx = key / 3;
    size_t k = idx % (g->nz + 1);
    size_t row = idx / (g->nz + 1);
    size_t j = row % (g->ny + 1);
    size_t i = row / (g->ny + 1);

    return i >= i0 && i < i1 &&
           j >= j0 && j < j1 &&
           k >= k0 && k < k1;
}

int materials_assign_box(
    materials_t *mat,
    const fdtd_grid_t *g,
    size_t i0,
    size_t i1,
    size_t j0,
    size_t j1,
    size_t k0,
    size_t k1,
    material_id_t id) {
    if (!box_valid(mat, g, i0, i1, j0, j1, k0, k1) ||
        (size_t) id >= mat->table_count) {
        return -1;
    }

    /* Clip periodic endpoints instead of wrapping them onto zero. */
    if ((g->periodic_axes & 1u) && i1 > g->nx) i1 = g->nx;
    if ((g->periodic_axes & 2u) && j1 > g->ny) j1 = g->ny;
    if ((g->periodic_axes & 4u) && k1 > g->nz) k1 = g->nz;

    for (size_t i = i0; i < i1; i++) {
        for (size_t j = j0; j < j1; j++) {
            for (size_t k = k0; k < k1; k++) {
                mat->ids[fdtd_index(g, i, j, k)] = id;
            }
        }
    }

    size_t kept = 0;

    for (size_t n = 0; n < mat->override_count; n++) {
        if (!key_inside_box(
            g, mat->overrides[n].key,
            i0, i1, j0, j1, k0, k1)) {
            mat->overrides[kept++] = mat->overrides[n];
        }
    }

    mat->override_count = kept;
    kept = 0;

    for (size_t n = 0; n < mat->tensor_override_count; n++) {
        if (!key_inside_box(
            g, mat->tensor_overrides[n].key,
            i0, i1, j0, j1, k0, k1)) {
            mat->tensor_overrides[kept++] =
                    mat->tensor_overrides[n];
        }
    }

    mat->tensor_override_count = kept;
    return 0;
}

int materials_set_box(
    materials_t *mat,
    const fdtd_grid_t *g,
    size_t i0,
    size_t i1,
    size_t j0,
    size_t j1,
    size_t k0,
    size_t k1,
    double eps_r,
    double sigma) {
    if (!box_valid(mat, g, i0, i1, j0, j1, k0, k1)) {
        return -1;
    }

    material_id_t id;

    if (materials_register(mat, eps_r, sigma, &id) != 0) {
        return -1;
    }

    return materials_assign_box(
        mat, g, i0, i1, j0, j1, k0, k1, id);
}

int materials_set_component(
    materials_t *mat,
    const fdtd_grid_t *g,
    field_component_t component,
    size_t i,
    size_t j,
    size_t k,
    double eps_r,
    double sigma) {
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

    if (reserve_override(mat) != 0) {
        return -1;
    }

    size_t trailing = mat->override_count - pos;

    if (trailing != 0) {
        memmove(
            mat->overrides + pos + 1,
            mat->overrides + pos,
            trailing * sizeof(*mat->overrides));
    }

    mat->overrides[pos] = (override_t){
        .key = key,
        .coeff = coeff
    };

    mat->override_count++;
    return 0;
}

int materials_clear_component(
    materials_t *mat,
    const fdtd_grid_t *g,
    field_component_t component,
    size_t i,
    size_t j,
    size_t k) {
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
                trailing * sizeof(*mat->overrides));
        }

        mat->override_count--;
    }

    return 0;
}

int materials_get_coefficients(
    const materials_t *mat,
    const fdtd_grid_t *g,
    field_component_t component,
    size_t i,
    size_t j,
    size_t k,
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

/*
 * Tensor entries take precedence over scalar entries during updates.
 * Scalar setters/clearers do not remove tensor entries.
 *
 * Component-local tensor convention:
 *   Ex: (m11, m12, m22) = (minv_xx, minv_xy, minv_yy)
 *   Ey: (m11, m12, m22) = (minv_yy, minv_xy, minv_xx)
 *
 * Current physical scope:
 *   - lossless, z-invariant TM;
 *   - tensor stencils outside CPML slabs.
 *
 * These physical restrictions are caller obligations. This setter
 * validates coordinates, interpolation margins, and pointwise SPD.
 */
int materials_set_component_tensor(
    materials_t *mat,
    const fdtd_grid_t *g,
    field_component_t component,
    size_t i,
    size_t j,
    size_t k,
    double m11,
    double m12,
    double m22) {
    int c = component_number(component);

    if (c != 0 && c != 1) {
        return -1;
    }

    if (!grid_matches(mat, g) ||
        !component_valid(g, c, i, j, k) ||
        !tensor_margin_ok(g, i, j) ||
        !tensor_spd(m11, m12, m22)) {
        return -1;
    }

    double scale = mat->dt / (EPS0 * mat->dx);

    if (!isfinite(scale) || scale <= 0.0) {
        return -1;
    }

    material_tensor_t coeff = {
        .m11 = m11,
        .m12 = m12,
        .m22 = m22,
        .scale = scale
    };

    size_t key = fdtd_index(g, i, j, k) * 3 + (size_t) c;
    size_t pos = tensor_lower_bound(mat, key);

    if (pos < mat->tensor_override_count &&
        mat->tensor_overrides[pos].key == key) {
        mat->tensor_overrides[pos].coeff = coeff;
        return 0;
    }

    if (tensor_reserve(mat) != 0) {
        return -1;
    }

    size_t trailing = mat->tensor_override_count - pos;

    if (trailing != 0) {
        memmove(
            mat->tensor_overrides + pos + 1,
            mat->tensor_overrides + pos,
            trailing * sizeof(*mat->tensor_overrides));
    }

    mat->tensor_overrides[pos] = (tensor_override_t){
        .key = key,
        .coeff = coeff
    };

    mat->tensor_override_count++;
    return 0;
}

int materials_clear_component_tensor(
    materials_t *mat,
    const fdtd_grid_t *g,
    field_component_t component,
    size_t i,
    size_t j,
    size_t k) {
    int c = component_number(component);

    if (!grid_matches(mat, g) ||
        !component_valid(g, c, i, j, k)) {
        return -1;
    }

    size_t key = fdtd_index(g, i, j, k) * 3 + (size_t) c;
    size_t pos = tensor_lower_bound(mat, key);

    if (pos < mat->tensor_override_count &&
        mat->tensor_overrides[pos].key == key) {
        size_t trailing = mat->tensor_override_count - pos - 1;

        if (trailing != 0) {
            memmove(
                mat->tensor_overrides + pos,
                mat->tensor_overrides + pos + 1,
                trailing * sizeof(*mat->tensor_overrides));
        }

        mat->tensor_override_count--;
    }

    return 0;
}

int materials_get_tensor(
    const materials_t *mat,
    const fdtd_grid_t *g,
    field_component_t component,
    size_t i,
    size_t j,
    size_t k,
    material_tensor_t *out) {
    int c = component_number(component);

    if (!out ||
        !grid_matches(mat, g) ||
        !component_valid(g, c, i, j, k)) {
        return -1;
    }

    const material_tensor_t *t =
            tensor_lookup(mat, fdtd_index(g, i, j, k), c);

    if (!t) {
        return 1;
    }

    *out = *t;
    return 0;
}

/*
 * Zero is reached on backward-difference axes only when periodic.
 * Tensor margin validation similarly protects nonperiodic neighbors.
 */
static inline size_t material_previous(size_t q, size_t extent) {
    return q == 0 ? extent - 1 : q - 1;
}

static inline size_t material_next(size_t q, size_t extent) {
    return q + 1 == extent ? 0 : q + 1;
}

#define EX(i, j, k) g->ex[fdtd_index(g, (i), (j), (k))]
#define EY(i, j, k) g->ey[fdtd_index(g, (i), (j), (k))]
#define EZ(i, j, k) g->ez[fdtd_index(g, (i), (j), (k))]
#define HX(i, j, k) g->hx[fdtd_index(g, (i), (j), (k))]
#define HY(i, j, k) g->hy[fdtd_index(g, (i), (j), (k))]
#define HZ(i, j, k) g->hz[fdtd_index(g, (i), (j), (k))]

/*
 * Shared tensor update for ordinary and CPML callers.
 *
 * Diagonal entries remain m11. Cross terms use reciprocal edge pairs
 * and normalized correlations derived from m11, m12, and m22.
 *
 * CPML callers use this only where the whole tensor stencil lies
 * in the unstretched interior. No anisotropic CPML is implemented.
 */
/*
 * Reciprocal tensor coupling.
 *
 * Cross-couplings exist only between two tensor-enabled edges.
 * Each Ex/Ey pair uses the same coefficient in both directions.
 * The diagonal entries remain component-local m11 values.
 */
static double tensor_pair_weight(
    const material_tensor_t *ex,
    const material_tensor_t *ey)
{
    double rx = (ex->m12 / sqrt(ex->m11)) / sqrt(ex->m22);
    double ry = (ey->m12 / sqrt(ey->m11)) / sqrt(ey->m22);
    double rho = 0.5 * rx + 0.5 * ry;

    return (0.25 * rho * sqrt(ex->m11)) * sqrt(ey->m11);
}

static inline void update_tensor_edge(
    fdtd_grid_t *g,
    const materials_t *mat,
    const material_tensor_t *t,
    int component,
    size_t i,
    size_t j,
    size_t k)
{
    double value;

    if (component == 0) {
        size_t ip = material_next(i, g->nx);
        size_t jm = material_previous(j, g->ny);
        size_t ni[4] = { i, ip, i, ip };
        size_t nj[4] = { j, j, jm, jm };

        value = t->m11 * (HZ(i, j, k) - HZ(i, jm, k));

        for (size_t q = 0; q < 4; q++) {
            size_t u = ni[q];
            size_t v = nj[q];
            const material_tensor_t *other = tensor_lookup(
                mat, fdtd_index(g, u, v, k), 1);

            if (!other) continue;

            size_t um = material_previous(u, g->nx);
            double curl_y = -(HZ(u, v, k) - HZ(um, v, k));
            value += tensor_pair_weight(t, other) * curl_y;
        }

        EX(i, j, k) += t->scale * value;
    } else {
        size_t im = material_previous(i, g->nx);
        size_t jp = material_next(j, g->ny);
        size_t ni[4] = { i, im, i, im };
        size_t nj[4] = { j, j, jp, jp };

        value = -t->m11 * (HZ(i, j, k) - HZ(im, j, k));

        for (size_t q = 0; q < 4; q++) {
            size_t u = ni[q];
            size_t v = nj[q];
            const material_tensor_t *other = tensor_lookup(
                mat, fdtd_index(g, u, v, k), 0);

            if (!other) continue;

            size_t vm = material_previous(v, g->ny);
            double curl_x = HZ(u, v, k) - HZ(u, vm, k);
            value += tensor_pair_weight(other, t) * curl_x;
        }

        EY(i, j, k) += t->scale * value;
    }
}

/*
 * A NULL pml selects ordinary updates.
 *
 * Keeping the component loops in one implementation ensures both
 * public update functions use the same tensor lookup and stencil.
 */
static void update_e_impl(
    fdtd_grid_t *g,
    cpml_t *pml,
    const materials_t *mat) {
    size_t i0 = (g->periodic_axes & 1u) ? 0u : 1u;
    size_t j0 = (g->periodic_axes & 2u) ? 0u : 1u;
    size_t k0 = (g->periodic_axes & 4u) ? 0u : 1u;

    /* Ex. */
    for (size_t i = 0; i < g->nx; i++) {
        for (size_t j = j0; j < g->ny; j++) {
            size_t jm = material_previous(j, g->ny);

            for (size_t k = k0; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);

                const material_tensor_t *t =
                        tensor_lookup(mat, idx, 0);

                if (t) {
                    update_tensor_edge(g, mat, t, 0, i, j, k);
                    continue;
                }

                size_t km = material_previous(k, g->nz);
                material_coeff_t coeff = lookup(mat, idx, 0);

                double termA =
                        HZ(i, j, k) - HZ(i, jm, k);

                double termB =
                        HY(i, j, k) - HY(i, j, km);

                if (pml) {
                    termA = cpml_term(
                        termA,
                        pml->kappa_inv_m_y[j],
                        pml->b_m_y[j],
                        pml->a_m_y[j],
                        g->dx,
                        &pml->psi_ex_y[idx]);

                    termB = cpml_term(
                        termB,
                        pml->kappa_inv_m_z[k],
                        pml->b_m_z[k],
                        pml->a_m_z[k],
                        g->dx,
                        &pml->psi_ex_z[idx]);
                }

                EX(i, j, k) =
                        coeff.Ca * EX(i, j, k) +
                        coeff.Cb * (termA - termB);
            }
        }
    }

    /* Ey. */
    for (size_t i = i0; i < g->nx; i++) {
        size_t im = material_previous(i, g->nx);

        for (size_t j = 0; j < g->ny; j++) {
            for (size_t k = k0; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);

                const material_tensor_t *t =
                        tensor_lookup(mat, idx, 1);

                if (t) {
                    update_tensor_edge(g, mat, t, 1, i, j, k);
                    continue;
                }

                size_t km = material_previous(k, g->nz);
                material_coeff_t coeff = lookup(mat, idx, 1);

                double termA =
                        HX(i, j, k) - HX(i, j, km);

                double termB =
                        HZ(i, j, k) - HZ(im, j, k);

                if (pml) {
                    termA = cpml_term(
                        termA,
                        pml->kappa_inv_m_z[k],
                        pml->b_m_z[k],
                        pml->a_m_z[k],
                        g->dx,
                        &pml->psi_ey_z[idx]);

                    termB = cpml_term(
                        termB,
                        pml->kappa_inv_m_x[i],
                        pml->b_m_x[i],
                        pml->a_m_x[i],
                        g->dx,
                        &pml->psi_ey_x[idx]);
                }

                EY(i, j, k) =
                        coeff.Ca * EY(i, j, k) +
                        coeff.Cb * (termA - termB);
            }
        }
    }

    /* Ez remains scalar. */
    for (size_t i = i0; i < g->nx; i++) {
        size_t im = material_previous(i, g->nx);

        for (size_t j = j0; j < g->ny; j++) {
            size_t jm = material_previous(j, g->ny);

            for (size_t k = 0; k < g->nz; k++) {
                size_t idx = fdtd_index(g, i, j, k);
                material_coeff_t coeff = lookup(mat, idx, 2);

                double termA =
                        HY(i, j, k) - HY(im, j, k);

                double termB =
                        HX(i, j, k) - HX(i, jm, k);

                if (pml) {
                    termA = cpml_term(
                        termA,
                        pml->kappa_inv_m_x[i],
                        pml->b_m_x[i],
                        pml->a_m_x[i],
                        g->dx,
                        &pml->psi_ez_x[idx]);

                    termB = cpml_term(
                        termB,
                        pml->kappa_inv_m_y[j],
                        pml->b_m_y[j],
                        pml->a_m_y[j],
                        g->dx,
                        &pml->psi_ez_y[idx]);
                }

                EZ(i, j, k) =
                        coeff.Ca * EZ(i, j, k) +
                        coeff.Cb * (termA - termB);
            }
        }
    }
}

void fdtd_update_e_materials(
    fdtd_grid_t *g,
    const materials_t *mat) {
    if (!grid_matches(mat, g)) {
        return;
    }

    update_e_impl(g, NULL, mat);
}

void fdtd_update_e_cpml_materials(
    fdtd_grid_t *g,
    cpml_t *pml,
    const materials_t *mat) {
    if (!grid_matches(mat, g) ||
        !pml || pml->dx != g->dx) {
        return;
    }

    update_e_impl(g, pml, mat);
}

#undef EX
#undef EY
#undef EZ
#undef HX
#undef HY
#undef HZ
