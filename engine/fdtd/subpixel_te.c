#include "subpixel_te.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

#define SPTE_CAPACITY 16

typedef struct {
    double x;
    double y;
} spte_point;

static int valid_rectangle(
    const subpixel_te_rectangle_t *rectangle)
{
    return rectangle &&
        isfinite(rectangle->center_x) &&
        isfinite(rectangle->center_y) &&
        isfinite(rectangle->size_u) &&
        rectangle->size_u > 0.0 &&
        isfinite(rectangle->size_v) &&
        rectangle->size_v > 0.0 &&
        isfinite(rectangle->angle_rad);
}

int subpixel_te_rectangle_edge_fraction(
    const subpixel_te_rectangle_t *rectangle,
    double sample_x,
    double sample_y,
    int edge_axis,
    double length,
    double *fraction_out
) {
    double rel_x, rel_y, u, v, t_lo, t_hi;
    double du, dv, a_lo, a_hi, b_lo, b_hi, lo, hi;

    if (!rectangle || !fraction_out) return SUBPIXEL_TE_INVALID_ARGUMENT;
    if (!isfinite(sample_x) || !isfinite(sample_y)) return SUBPIXEL_TE_INVALID_ARGUMENT;
    if (!(edge_axis == 1 || edge_axis == 2)) return SUBPIXEL_TE_INVALID_ARGUMENT;
    if (!isfinite(length) || !(length > 0.0)) return SUBPIXEL_TE_INVALID_ARGUMENT;
    if (!isfinite(rectangle->center_x) || !isfinite(rectangle->center_y) ||
        !isfinite(rectangle->size_u) || !(rectangle->size_u >= 0.0) ||
        !isfinite(rectangle->size_v) || !(rectangle->size_v >= 0.0) ||
        !isfinite(rectangle->angle_rad)) return SUBPIXEL_TE_INVALID_ARGUMENT;

    rel_x = sample_x - rectangle->center_x;
    rel_y = sample_y - rectangle->center_y;
    if (!isfinite(rel_x) || !isfinite(rel_y)) return SUBPIXEL_TE_NUMERIC_ERROR;

    {
        double ca = cos(rectangle->angle_rad), sa = sin(rectangle->angle_rad);
        double dxu = (edge_axis == 1) ? ca : -sa;   /* du/dt or dv/dt sign-checked below */
        u  = rel_x * ca + rel_y * sa;
        v  = rel_y * ca - rel_x * sa;
        du = (edge_axis == 1) ? ca  : -sa;          /* u(t) = u + du*t, t in [0,length] */
        dv = (edge_axis == 1) ? sa  :  ca;          /* v(t) = v + dv*t                  */
        (void)dxu;
    }

    /* t interval where |u + du t| <= size_u/2 */
    if (!rectangle->size_u) return SUBPIXEL_TE_NUMERIC_ERROR;
    if (fabs(du) < 1e-300) {
        if (fabs(u) > 0.5 * rectangle->size_u * (1.0 + 1e-12)) { *fraction_out = 0.0; return SUBPIXEL_TE_OK; }
        a_lo = 0.0; a_hi = length;
    } else {
        t_lo = (-0.5 * rectangle->size_u - u) / du;
        t_hi = ( 0.5 * rectangle->size_u - u) / du;
        a_lo = fmin(t_lo, t_hi); a_hi = fmax(t_lo, t_hi);
    }

    /* same clip for v */
    if (fabs(dv) < 1e-300) {
        if (fabs(v) > 0.5 * rectangle->size_v * (1.0 + 1e-12)) { *fraction_out = 0.0; return SUBPIXEL_TE_OK; }
        b_lo = 0.0; b_hi = length;
    } else {
        t_lo = (-0.5 * rectangle->size_v - v) / dv;
        t_hi = ( 0.5 * rectangle->size_v - v) / dv;
        b_lo = fmin(t_lo, t_hi); b_hi = fmax(t_lo, t_hi);
    }

    lo = fmax(fmax(a_lo, b_lo), 0.0);
    hi = fmin(fmin(a_hi, b_hi), length);
    *fraction_out = hi > lo ? (hi - lo) / length : 0.0;
    return SUBPIXEL_TE_OK;
}

int subpixel_te_mix_epsilon_harmonic(
    double fraction,
    double eps_background,
    double eps_inside,
    double *epsilon_r_out
) {
    double inv;
    if (!epsilon_r_out) return SUBPIXEL_TE_INVALID_ARGUMENT;
    if (!isfinite(fraction) || fraction < 0.0 || fraction > 1.0) return SUBPIXEL_TE_INVALID_ARGUMENT;
    if (!isfinite(eps_background) || eps_background < 1.0) return SUBPIXEL_TE_INVALID_ARGUMENT;
    if (!isfinite(eps_inside) || eps_inside < 1.0) return SUBPIXEL_TE_INVALID_ARGUMENT;

    inv = (1.0 - fraction) / eps_background + fraction / eps_inside;
    if (!isfinite(inv) || !(inv > 0.0)) return SUBPIXEL_TE_NUMERIC_ERROR;
    *epsilon_r_out = 1.0 / inv;
    return SUBPIXEL_TE_OK;
}

int subpixel_te_mix_epsilon_tensor(
    double fraction,
    double eps_background,
    double eps_inside,
    double normal_rad,
    double *minv_xx, double *minv_xy, double *minv_yy
) {
    double lam_n, eps_t, lam_t, ca, sa;
    if (!minv_xx || !minv_xy || !minv_yy) return SUBPIXEL_TE_INVALID_ARGUMENT;
    if (!isfinite(fraction) || fraction < 0.0 || fraction > 1.0) return SUBPIXEL_TE_INVALID_ARGUMENT;
    if (!isfinite(eps_background) || eps_background < 1.0) return SUBPIXEL_TE_INVALID_ARGUMENT;
    if (!isfinite(eps_inside) || eps_inside < 1.0) return SUBPIXEL_TE_INVALID_ARGUMENT;
    if (!isfinite(normal_rad)) return SUBPIXEL_TE_INVALID_ARGUMENT;

    lam_n = (1.0 - fraction) / eps_background + fraction / eps_inside;
    eps_t = (1.0 - fraction) * eps_background + fraction * eps_inside;
    if (!isfinite(lam_n) || !(lam_n > 0.0) || !isfinite(eps_t) || !(eps_t > 0.0))
        return SUBPIXEL_TE_NUMERIC_ERROR;
    lam_t = 1.0 / eps_t;

    ca = cos(normal_rad); sa = sin(normal_rad);

    *minv_xx = lam_n * ca * ca + lam_t * sa * sa;
    *minv_yy = lam_n * sa * sa + lam_t * ca * ca;
    *minv_xy = (lam_n - lam_t) * ca * sa;
    return SUBPIXEL_TE_OK;
}

static int append_point(
    spte_point *output,
    size_t *count,
    spte_point point)
{
    if (*count >= SPTE_CAPACITY ||
        !isfinite(point.x) ||
        !isfinite(point.y)) {
        return SUBPIXEL_TE_NUMERIC_ERROR;
    }

    output[(*count)++] = point;

    return SUBPIXEL_TE_OK;
}

static int clip_plane(
    const spte_point *input,
    size_t count,
    double a,
    double b,
    double d,
    spte_point *output,
    size_t *output_count)
{
    *output_count = 0;

    if (count == 0) {
        return SUBPIXEL_TE_OK;
    }

    if (!isfinite(a) ||
        !isfinite(b) ||
        !isfinite(d)) {
        return SUBPIXEL_TE_NUMERIC_ERROR;
    }

    double scale = fmax(
        fabs(a),
        fmax(fabs(b), fabs(d))
    );

    if (scale == 0.0) {
        for (size_t i = 0; i < count; i++) {
            output[i] = input[i];
        }

        *output_count = count;

        return SUBPIXEL_TE_OK;
    }

    a /= scale;
    b /= scale;
    d /= scale;

    spte_point previous = input[count - 1];

    double previous_distance =
        fma(a, previous.x, b * previous.y) - d;

    int previous_inside =
        previous_distance <= 0.0;

    for (size_t i = 0; i < count; i++) {
        spte_point current = input[i];

        double distance =
            fma(a, current.x, b * current.y) - d;

        int inside = distance <= 0.0;

        if (inside != previous_inside) {
            double denominator =
                previous_distance - distance;

            if (!isfinite(denominator) ||
                denominator == 0.0) {
                return SUBPIXEL_TE_NUMERIC_ERROR;
            }

            double t =
                previous_distance / denominator;

            if (!isfinite(t)) {
                return SUBPIXEL_TE_NUMERIC_ERROR;
            }

            t = fmax(0.0, fmin(1.0, t));

            spte_point crossing = {
                fma(
                    t,
                    current.x - previous.x,
                    previous.x
                ),
                fma(
                    t,
                    current.y - previous.y,
                    previous.y
                )
            };

            int status = append_point(
                output,
                output_count,
                crossing
            );

            if (status != SUBPIXEL_TE_OK) {
                return status;
            }
        }

        if (inside) {
            int status = append_point(
                output,
                output_count,
                current
            );

            if (status != SUBPIXEL_TE_OK) {
                return status;
            }
        }

        previous = current;
        previous_distance = distance;
        previous_inside = inside;
    }

    return SUBPIXEL_TE_OK;
}

int subpixel_te_rectangle_fraction(
    const subpixel_te_rectangle_t *rectangle,
    double sample_x,
    double sample_y,
    double dx,
    double dy,
    double *fraction_out)
{
    if (!fraction_out ||
        !valid_rectangle(rectangle) ||
        !isfinite(sample_x) ||
        !isfinite(sample_y) ||
        !isfinite(dx) ||
        dx <= 0.0 ||
        !isfinite(dy) ||
        dy <= 0.0) {
        return SUBPIXEL_TE_INVALID_ARGUMENT;
    }

    double half_u =
        0.5 * rectangle->size_u;

    double half_v =
        0.5 * rectangle->size_v;

    double relative_x =
        sample_x - rectangle->center_x;

    double relative_y =
        sample_y - rectangle->center_y;

    if (!(half_u > 0.0) ||
        !(half_v > 0.0) ||
        !isfinite(relative_x) ||
        !isfinite(relative_y)) {
        return SUBPIXEL_TE_NUMERIC_ERROR;
    }

    double c = cos(rectangle->angle_rad);
    double s = sin(rectangle->angle_rad);

    double center_u = fma(c, relative_x, s * relative_y);

    double center_v = fma(-s, relative_x, c * relative_y);

    double planes[4][3] = {
        {
            c * dx,
            s * dy,
            half_u - center_u
        },
        {
            -c * dx,
            -s * dy,
            half_u + center_u
        },
        {
            -s * dx,
            c * dy,
            half_v - center_v
        },
        {
            s * dx,
            -c * dy,
            half_v + center_v
        }
    };

    spte_point polygons[2][SPTE_CAPACITY] = {
        {
            {-0.5, -0.5},
            { 0.5, -0.5},
            { 0.5,  0.5},
            {-0.5,  0.5}
        },
        {
            {0.0, 0.0}
        }
    };

    size_t count = 4;
    int active = 0;

    for (size_t p = 0; p < 4; p++) {
        size_t next_count = 0;

        int status = clip_plane(
            polygons[active],
            count,
            planes[p][0],
            planes[p][1],
            planes[p][2],
            polygons[1 - active],
            &next_count
        );

        if (status != SUBPIXEL_TE_OK) {
            return status;
        }

        count = next_count;
        active = 1 - active;
    }

    double twice_area = 0.0;

    if (count >= 3) {
        const spte_point *polygon =
            polygons[active];

        for (size_t i = 1; i + 1 < count; i++) {
            double ax =
                polygon[i].x - polygon[0].x;

            double ay =
                polygon[i].y - polygon[0].y;

            double bx =
                polygon[i + 1].x - polygon[0].x;

            double by =
                polygon[i + 1].y - polygon[0].y;

            twice_area +=
                fma(ax, by, -ay * bx);
        }
    }

    double fraction =
        0.5 * fabs(twice_area);

    if (!isfinite(fraction) ||
        fraction > 1.0 + 128.0 * DBL_EPSILON) {
        return SUBPIXEL_TE_NUMERIC_ERROR;
    }

    *fraction_out =
        fmax(0.0, fmin(1.0, fraction));

    return SUBPIXEL_TE_OK;
}

int subpixel_te_mix_epsilon(
    double fraction,
    double eps_background,
    double eps_inside,
    double *epsilon_r_out)
{
    if (!epsilon_r_out ||
        !isfinite(fraction) ||
        fraction < 0.0 ||
        fraction > 1.0 ||
        !isfinite(eps_background) ||
        eps_background < 1.0 ||
        !isfinite(eps_inside) ||
        eps_inside < 1.0) {
        return SUBPIXEL_TE_INVALID_ARGUMENT;
    }

    double epsilon;

    if (fraction == 0.0) {
        epsilon = eps_background;
    } else if (fraction == 1.0) {
        epsilon = eps_inside;
    } else if (eps_inside >= eps_background) {
        epsilon = fma(
            fraction,
            eps_inside - eps_background,
            eps_background
        );
    } else {
        epsilon = fma(
            1.0 - fraction,
            eps_background - eps_inside,
            eps_inside
        );
    }

    if (!isfinite(epsilon) ||
        epsilon < 1.0) {
        return SUBPIXEL_TE_NUMERIC_ERROR;
    }

    *epsilon_r_out = epsilon;

    return SUBPIXEL_TE_OK;
}

int subpixel_te_rectangle_epsilon(
    const subpixel_te_rectangle_t *rectangle,
    double sample_x,
    double sample_y,
    double dx,
    double dy,
    double eps_background,
    double eps_inside,
    double *fraction_out,
    double *epsilon_r_out)
{
    if (!fraction_out ||
        !epsilon_r_out ||
        fraction_out == epsilon_r_out) {
        return SUBPIXEL_TE_INVALID_ARGUMENT;
    }

    double fraction;
    double epsilon;

    int status = subpixel_te_rectangle_fraction(
        rectangle,
        sample_x,
        sample_y,
        dx,
        dy,
        &fraction
    );

    if (status != SUBPIXEL_TE_OK) {
        return status;
    }

    status = subpixel_te_mix_epsilon(
        fraction,
        eps_background,
        eps_inside,
        &epsilon
    );

    if (status != SUBPIXEL_TE_OK) {
        return status;
    }

    *fraction_out = fraction;
    *epsilon_r_out = epsilon;

    return SUBPIXEL_TE_OK;
}