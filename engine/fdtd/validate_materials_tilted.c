#include "extended_test.h"
#include "subpixel_te.h"

#define TILTED_PROBES 4

typedef enum {
    SMOOTHING_OFF = 0,
    SMOOTHING_TE,
    SMOOTHING_TENSOR
} smoothing_mode_t;

static const double probe_x[TILTED_PROBES] = {
    0.030, 0.030, 0.030, 0.062
};

static const double probe_y[TILTED_PROBES] = {
    0.032, 0.040, 0.048, 0.040
};

static const char *smoothing_label(smoothing_mode_t mode)
{
    switch (mode) {
        case SMOOTHING_TE:
            return "te";
        case SMOOTHING_TENSOR:
            return "tensor";
        default:
            return "off";
    }
}

static int inside_slab(
    double x,
    double y,
    double angle,
    double normal_shift)
{
    double u =
        (x - 0.045) * cos(angle) +
        (y - 0.040) * sin(angle) -
        normal_shift;

    double v =
        -(x - 0.045) * sin(angle) +
        (y - 0.040) * cos(angle);

    return fabs(u) < 0.003 && fabs(v) < 0.010;
}

/*
 * Clip one coordinate of a segment against [-half_size, half_size].
 * The segment parameter is t in [0, 1].
 */
static int clip_segment_axis(
    double start,
    double delta,
    double half_size,
    double *t0,
    double *t1)
{
    if (delta == 0.0) {
        return fabs(start) <= half_size;
    }

    double a = (-half_size - start) / delta;
    double b = ( half_size - start) / delta;

    if (a > b) {
        double tmp = a;
        a = b;
        b = tmp;
    }

    if (a > *t0) {
        *t0 = a;
    }

    if (b < *t1) {
        *t1 = b;
    }

    return *t1 >= *t0;
}

static int rectangle_edge_coverage(
    const subpixel_te_rectangle_t *rectangle,
    double x,
    double y,
    int axis,
    double length,
    double *fraction_out)
{
    if (!rectangle || !fraction_out ||
        (axis != 0 && axis != 1) ||
        !isfinite(length) || length <= 0.0) {
        return -1;
    }

    double c = cos(rectangle->angle_rad);
    double s = sin(rectangle->angle_rad);

    double rx = x - rectangle->center_x;
    double ry = y - rectangle->center_y;

    double u0 =  c * rx + s * ry;
    double v0 = -s * rx + c * ry;

    double du = axis == 0 ? c * length : s * length;
    double dv = axis == 0 ? -s * length : c * length;

    double t0 = 0.0;
    double t1 = 1.0;

    if (!clip_segment_axis(
            u0, du, 0.5 * rectangle->size_u, &t0, &t1) ||
        !clip_segment_axis(
            v0, dv, 0.5 * rectangle->size_v, &t0, &t1)) {
        *fraction_out = 0.0;
        return 0;
    }

    double fraction = t1 - t0;

    if (!isfinite(fraction)) {
        return -1;
    }

    if (fraction < 0.0) {
        fraction = 0.0;
    } else if (fraction > 1.0) {
        fraction = 1.0;
    }

    *fraction_out = fraction;
    return 0;
}

static double rectangle_normal_angle(
    const subpixel_te_rectangle_t *rectangle,
    double x,
    double y)
{
    double c = cos(rectangle->angle_rad);
    double s = sin(rectangle->angle_rad);

    double rx = x - rectangle->center_x;
    double ry = y - rectangle->center_y;

    double u =  c * rx + s * ry;
    double v = -s * rx + c * ry;

    double distance_u =
        fabs(fabs(u) - 0.5 * rectangle->size_u);

    double distance_v =
        fabs(fabs(v) - 0.5 * rectangle->size_v);

    return distance_u <= distance_v
        ? rectangle->angle_rad
        : rectangle->angle_rad + 0.5 * XT_PI;
}

static int assign_tensor_edge(
    xt_sim *s,
    const subpixel_te_rectangle_t *rectangle,
    size_t i,
    size_t j,
    int axis)
{
    double dx = s->g->dx;
    double x = i * dx;
    double y = j * dx;
    double fraction = 0.0;

    if (rectangle_edge_coverage(
            rectangle, x, y, axis, dx, &fraction) != 0) {
        fprintf(
            stderr,
            "Tensor edge coverage failed at (%zu, %zu), axis %d\n",
            i, j, axis);
        return -1;
    }

    if (fraction <= 1e-12 || fraction >= 1.0 - 1e-12) {
        return 0;
    }

    double midpoint_x = x + (axis == 0 ? 0.5 * dx : 0.0);
    double midpoint_y = y + (axis == 1 ? 0.5 * dx : 0.0);

    double normal_angle =
        rectangle_normal_angle(rectangle, midpoint_x, midpoint_y);

    double mxx = 0.0;
    double mxy = 0.0;
    double myy = 0.0;

    int status = subpixel_te_mix_epsilon_tensor(
        fraction,
        1.0,
        4.0,
        normal_angle,
        &mxx,
        &mxy,
        &myy);

    if (status != SUBPIXEL_TE_OK) {
        fprintf(
            stderr,
            "Tensor mixing failed at (%zu, %zu), axis %d, status %d\n",
            i, j, axis, status);
        return -1;
    }

    field_component_t component = axis == 0 ? FIELD_EX : FIELD_EY;

    /*
     * Current materials API/update convention:
     * m11 is the aligned coefficient at BOTH Ex and Ey edges.
     * Thus Ey stores yy first and xx second.
     */
    double aligned = axis == 0 ? mxx : myy;
    double transverse = axis == 0 ? myy : mxx;

    if (materials_set_component_tensor(
            s->m,
            s->g,
            component,
            i, j, 0,
            aligned,
            mxy,
            transverse) != 0) {
        fprintf(
            stderr,
            "Tensor assignment failed at (%zu, %zu), axis %d\n",
            i, j, axis);
        return -1;
    }

    return 0;
}

static int build_geometry(
    xt_sim *s,
    double angle,
    double normal_shift,
    smoothing_mode_t smoothing)
{
    material_id_t dielectric;

    subpixel_te_rectangle_t rectangle = {
        .center_x = 0.045 + normal_shift * cos(angle),
        .center_y = 0.040 + normal_shift * sin(angle),
        .size_u = 0.006,
        .size_v = 0.020,
        .angle_rad = angle
    };

    if (materials_register(s->m, 4.0, 0.0, &dielectric) != 0) {
        fprintf(stderr, "Could not register dielectric material\n");
        return -1;
    }

    for (size_t i = 0; i < s->g->nx; i++) {
        for (size_t j = 0; j < s->g->ny; j++) {
            double x = i * s->g->dx;
            double y = j * s->g->dx;

            int base_inside =
                inside_slab(x, y, angle, normal_shift);

            int ex_inside =
                inside_slab(
                    x + 0.5 * s->g->dx,
                    y,
                    angle,
                    normal_shift);

            int ey_inside =
                inside_slab(
                    x,
                    y + 0.5 * s->g->dx,
                    angle,
                    normal_shift);

            if (base_inside) {
                if (materials_assign_box(
                        s->m, s->g,
                        i, i + 1,
                        j, j + 1,
                        0, 2,
                        dielectric) != 0) {
                    fprintf(
                        stderr,
                        "Material assignment failed at (%zu, %zu)\n",
                        i, j);
                    return -1;
                }
            }

            if (ex_inside != base_inside) {
                if (materials_set_component(
                        s->m, s->g, FIELD_EX,
                        i, j, 0,
                        ex_inside ? 4.0 : 1.0,
                        0.0) != 0) {
                    fprintf(
                        stderr,
                        "Ex material assignment failed at (%zu, %zu)\n",
                        i, j);
                    return -1;
                }
            }

            if (ey_inside != base_inside) {
                if (materials_set_component(
                        s->m, s->g, FIELD_EY,
                        i, j, 0,
                        ey_inside ? 4.0 : 1.0,
                        0.0) != 0) {
                    fprintf(
                        stderr,
                        "Ey material assignment failed at (%zu, %zu)\n",
                        i, j);
                    return -1;
                }
            }

            if (smoothing == SMOOTHING_TE) {
                double fraction;
                double epsilon;
                double base_epsilon = base_inside ? 4.0 : 1.0;

                int status = subpixel_te_rectangle_epsilon(
                    &rectangle,
                    x, y,
                    s->g->dx, s->g->dx,
                    1.0, 4.0,
                    &fraction,
                    &epsilon);

                if (status != SUBPIXEL_TE_OK) {
                    fprintf(
                        stderr,
                        "TE smoothing failed at (%zu, %zu), status %d\n",
                        i, j, status);
                    return -1;
                }

                if (epsilon != base_epsilon) {
                    if (materials_set_component(
                            s->m, s->g, FIELD_EZ,
                            i, j, 0,
                            epsilon, 0.0) != 0) {
                        fprintf(
                            stderr,
                            "Ez smoothing assignment failed at (%zu, %zu)\n",
                            i, j);
                        return -1;
                    }
                }
            }
        }
    }

    if (smoothing == SMOOTHING_TENSOR) {
        for (size_t i = 0; i < s->g->nx; i++) {
            for (size_t j = 0; j < s->g->ny; j++) {
                if (assign_tensor_edge(
                        s, &rectangle, i, j, 0) != 0 ||
                    assign_tensor_edge(
                        s, &rectangle, i, j, 1) != 0) {
                    return -1;
                }
            }
        }
    }

    return 0;
}

static double sample_probe(
    fdtd_grid_t *g,
    int tm,
    double x,
    double y)
{
    double gx = x / g->dx;
    double gy = y / g->dx - (tm ? 0.5 : 0.0);

    size_t i = (size_t)floor(gx);
    size_t j = (size_t)floor(gy);

    double a = gx - i;
    double b = gy - j;
    double *field = xt_e(g, tm);

    return
        (1.0 - a) * (1.0 - b) *
            field[fdtd_index(g, i, j, 0)] +
        a * (1.0 - b) *
            field[fdtd_index(g, i + 1, j, 0)] +
        (1.0 - a) * b *
            field[fdtd_index(g, i, j + 1, 0)] +
        a * b *
            field[fdtd_index(g, i + 1, j + 1, 0)];
}

static int parse_number(const char *text, double *value)
{
    char *end = NULL;

    errno = 0;
    *value = strtod(text, &end);

    return errno != 0 ||
           end == text ||
           *end != '\0' ||
           !isfinite(*value);
}

int main(int argc, char **argv)
{
    int scale = RFSYN_VALIDATION_REFINEMENT;
    int tm = 0;
    smoothing_mode_t smoothing = SMOOTHING_OFF;

    int settle_cycles = 24;
    int window_cycles = 12;
    int max_windows = 8;

    double angle_deg = 30.0;
    double shift_mm = 0.5;
    const char *output_name = NULL;

    for (int q = 1; q < argc; q++) {
        if (strcmp(argv[q], "--help") == 0) {
            printf(
                "Usage: %s "
                "[--pol TE|TM] "
                "[--smoothing off|te|tensor] "
                "[--angle 0..60] "
                "[--shift-mm -1..1] "
                "[--scale 1|2|4|8] "
                "[--settle-cycles 8..500] "
                "[--window-cycles 4..100] "
                "[--max-windows 3..100] "
                "[--output filename]\n",
                argv[0]);
            return 0;
        }

        if (q + 1 >= argc) {
            fprintf(stderr, "Missing value for %s\n", argv[q]);
            return 2;
        }

        const char *key = argv[q++];
        const char *text = argv[q];

        if (strcmp(key, "--pol") == 0) {
            if (strcmp(text, "TE") == 0) {
                tm = 0;
            } else if (strcmp(text, "TM") == 0) {
                tm = 1;
            } else {
                fprintf(stderr, "--pol must be TE or TM\n");
                return 2;
            }
            continue;
        }

        if (strcmp(key, "--smoothing") == 0) {
            if (strcmp(text, "off") == 0) {
                smoothing = SMOOTHING_OFF;
            } else if (strcmp(text, "te") == 0) {
                smoothing = SMOOTHING_TE;
            } else if (strcmp(text, "tensor") == 0) {
                smoothing = SMOOTHING_TENSOR;
            } else {
                fprintf(
                    stderr,
                    "--smoothing must be off, te, or tensor\n");
                return 2;
            }
            continue;
        }

        if (strcmp(key, "--output") == 0) {
            if (*text == '\0') {
                fprintf(stderr, "--output cannot be empty\n");
                return 2;
            }
            output_name = text;
            continue;
        }

        double value;

        if (parse_number(text, &value) != 0) {
            fprintf(
                stderr,
                "Invalid numeric value: %s %s\n",
                key, text);
            return 2;
        }

        if (strcmp(key, "--angle") == 0 &&
            value >= 0.0 && value <= 60.0) {
            angle_deg = value;
        } else if (strcmp(key, "--shift-mm") == 0 &&
                   fabs(value) <= 1.0) {
            shift_mm = value;
        } else if (strcmp(key, "--scale") == 0 &&
                   (value == 1.0 || value == 2.0 ||
                    value == 4.0 || value == 8.0)) {
            scale = (int)value;
        } else if (strcmp(key, "--settle-cycles") == 0 &&
                   value >= 8.0 && value <= 500.0 &&
                   value == floor(value)) {
            settle_cycles = (int)value;
        } else if (strcmp(key, "--window-cycles") == 0 &&
                   value >= 4.0 && value <= 100.0 &&
                   value == floor(value)) {
            window_cycles = (int)value;
        } else if (strcmp(key, "--max-windows") == 0 &&
                   value >= 3.0 && value <= 100.0 &&
                   value == floor(value)) {
            max_windows = (int)value;
        } else {
            fprintf(stderr, "Invalid option: %s %s\n", key, text);
            return 2;
        }
    }

    if (scale != 1 && scale != 2 && scale != 4 && scale != 8) {
        fprintf(stderr, "Scale must be 1, 2, 4, or 8\n");
        return 2;
    }

    if (tm && smoothing == SMOOTHING_TE) {
        fprintf(stderr, "--smoothing te only supports --pol TE\n");
        return 2;
    }

    if (!tm && smoothing == SMOOTHING_TENSOR) {
        fprintf(stderr, "--smoothing tensor only supports --pol TM\n");
        return 2;
    }

    if (smoothing == SMOOTHING_TENSOR) {
        fprintf(
            stderr,
            "Experimental tensor mode: edge coverage with a "
            "nearest-face normal; accuracy is not yet validated.\n");
    }

    if (!output_name) {
        if (smoothing == SMOOTHING_TE) {
            output_name = "materials_tilted_te_samples.csv";
        } else if (smoothing == SMOOTHING_TENSOR) {
            output_name = "materials_tilted_tm_tensor_samples.csv";
        } else {
            output_name = "materials_tilted_samples.csv";
        }
    }

    const char *smoothing_name = smoothing_label(smoothing);

    xt_sim sims[2] = {{0}};
    FILE *csv = NULL;
    int result = 2;

    double dx = 0.001 / scale;
    double angle = angle_deg * XT_PI / 180.0;

    if (xt_create(
            &sims[0],
            80 * scale, 80 * scale, 1,
            dx, 4, 12 * scale) != 0) {
        fprintf(stderr, "Could not create vacuum reference simulation\n");
        goto cleanup;
    }

    if (xt_create(
            &sims[1],
            80 * scale, 80 * scale, 1,
            dx, 4, 12 * scale) != 0) {
        fprintf(stderr, "Could not create dielectric simulation\n");
        goto cleanup;
    }

    if (build_geometry(
            &sims[1],
            angle,
            shift_mm * 0.001,
            smoothing) != 0) {
        goto cleanup;
    }

    const double frequency = C0 / 0.024;
    double omega_dt = 2.0 * XT_PI * frequency * sims[0].g->dt;

    size_t period = (size_t)ceil(2.0 * XT_PI / omega_dt);
    size_t ramp_steps = 8 * period;
    size_t settle_steps = (size_t)settle_cycles * period;
    size_t window_steps = (size_t)window_cycles * period;
    size_t maximum_steps =
        settle_steps + (size_t)max_windows * window_steps;

    double complex previous[2][TILTED_PROBES] = {{0}};
    double complex current[2][TILTED_PROBES] = {{0}};

    double window_drift = INFINITY;
    double incident_norm = 0.0;
    int converged = 0;
    size_t final_step = 0;

    for (size_t n = 0; n < maximum_steps; n++) {
        for (size_t g = 0; g < 2; g++) {
            xt_step(&sims[g]);

            fdtd_grid_t *grid = sims[g].g;

            double envelope = n < ramp_steps
                ? pow(
                    sin(
                        0.5 * XT_PI *
                        (n + 1) / ramp_steps),
                    2.0)
                : 1.0;

            for (size_t j = 1; j < grid->ny; j++) {
                double y = (j + (tm ? 0.5 : 0.0)) * dx;

                if (y >= 0.016 && y <= 0.064) {
                    double source_value =
                        envelope *
                        cos(omega_dt * (n + 1)) *
                        exp(-pow((y - 0.040) / 0.0075, 2.0));

                    size_t idx =
                        fdtd_index(grid, 18 * scale, j, 0);

                    xt_e(grid, tm)[idx] += source_value;
                }
            }
        }

        if (n < settle_steps) {
            continue;
        }

        size_t sample_in_window =
            (n - settle_steps) % window_steps;

        size_t window_number =
            (n - settle_steps) / window_steps;

        double weight =
            4.0 / window_steps *
            pow(
                sin(
                    XT_PI *
                    (sample_in_window + 0.5) /
                    window_steps),
                2.0);

        double complex temporal_weight =
            weight * cexp(I * omega_dt * (n + 1));

        for (size_t g = 0; g < 2; g++) {
            for (size_t p = 0; p < TILTED_PROBES; p++) {
                double value = sample_probe(
                    sims[g].g,
                    tm,
                    probe_x[p],
                    probe_y[p]);

                if (!isfinite(value)) {
                    fprintf(
                        stderr,
                        "Nonfinite monitored field at step %zu\n",
                        n);
                    goto cleanup;
                }

                current[g][p] += value * temporal_weight;
            }
        }

        if (sample_in_window + 1 != window_steps) {
            continue;
        }

        double norm_squared = 0.0;
        double delta_squared = 0.0;

        for (size_t p = 0; p < TILTED_PROBES; p++) {
            norm_squared += pow(cabs(current[0][p]), 2.0);
        }

        if (!isfinite(norm_squared) || norm_squared <= 1e-20) {
            fprintf(
                stderr,
                "Incident reference is invalid or too small\n");
            goto cleanup;
        }

        for (size_t g = 0; g < 2; g++) {
            for (size_t p = 0; p < TILTED_PROBES; p++) {
                if (!xt_finite_complex(current[g][p])) {
                    fprintf(
                        stderr,
                        "Nonfinite accumulated phasor at step %zu\n",
                        n);
                    goto cleanup;
                }

                delta_squared += pow(
                    cabs(current[g][p] - previous[g][p]),
                    2.0);
            }
        }

        if (!isfinite(delta_squared)) {
            fprintf(
                stderr,
                "Nonfinite convergence diagnostic at step %zu\n",
                n);
            goto cleanup;
        }

        incident_norm = sqrt(norm_squared);
        window_drift = sqrt(delta_squared / norm_squared);
        final_step = n + 1;

        if (window_number >= 1 && window_drift <= 1e-6) {
            converged = 1;
            break;
        }

        memcpy(previous, current, sizeof(previous));

        if (window_number + 1 < (size_t)max_windows) {
            memset(current, 0, sizeof(current));
        }
    }

    csv = fopen(output_name, "w");

    if (!csv) {
        perror(output_name);
        goto cleanup;
    }

    if (fprintf(
            csv,
            "angle_deg,shift_mm,scale,polarization,dx_m,"
            "frequency_hz,steps,probe_id,x_m,y_m,"
            "incident_re,incident_im,scattered_re,scattered_im,"
            "window_drift,incident_norm,measurement_converged,"
            "smoothing\n") < 0) {
        goto cleanup;
    }

    for (size_t p = 0; p < TILTED_PROBES; p++) {
        double complex incident = current[0][p];
        double complex scattered = current[1][p] - incident;

        if (!xt_finite_complex(scattered)) {
            fprintf(
                stderr,
                "Nonfinite scattered phasor at probe %zu\n",
                p);
            goto cleanup;
        }

        if (fprintf(
                csv,
                "%.17g,%.17g,%d,%s,%.17g,%.17g,%zu,%zu,"
                "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                "%.17g,%.17g,%d,%s\n",
                angle_deg,
                shift_mm,
                scale,
                tm ? "TM" : "TE",
                dx,
                frequency,
                final_step,
                p,
                probe_x[p],
                probe_y[p],
                creal(incident),
                cimag(incident),
                creal(scattered),
                cimag(scattered),
                window_drift,
                incident_norm,
                converged,
                smoothing_name) < 0) {
            goto cleanup;
        }
    }

    result = converged ? 0 : 1;

    printf(
        "Tilted finite slab: %s smoothing=%s "
        "angle=%.1f shift=%.3f mm scale=%d drift=%.6g %s\n",
        tm ? "TM" : "TE",
        smoothing_name,
        angle_deg,
        shift_mm,
        scale,
        window_drift,
        converged
            ? "MEASUREMENT CONVERGED (accuracy not yet graded)"
            : "UNSETTLED");

cleanup:
    if (csv && fclose(csv) != 0) {
        perror("Closing tilted measurement output");
        result = 2;
    }

    xt_destroy(&sims[0]);
    xt_destroy(&sims[1]);

    if (result == 2) {
        fprintf(
            stderr,
            "Tilted measurement failed: setup, field, or file error\n");
    }

    return result;
}