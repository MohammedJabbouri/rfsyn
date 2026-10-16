#include "extended_test.h"

static double complex project_field(
    fdtd_grid_t *g,
    size_t i_plane,
    double ky,
    int tm)
{
    double complex sum = 0.0;
    double *field = xt_e(g, tm);

    for (size_t j = 0; j < g->ny; j++) {
        double y_position = j + (tm ? 0.5 : 0.0);
        size_t idx = fdtd_index(g, i_plane, j, 0);

        sum += field[idx] * cexp(-I * ky * y_position);
    }

    return sum / g->ny;
}

static int run_case(
    FILE *csv,
    double requested_angle,
    int tm,
    int full)
{
    size_t scale = RFSYN_VALIDATION_REFINEMENT;

    size_t ny = (full ? 64 : 32) * scale;
    size_t nx_test = (full ? 256 : 128) * scale;
    size_t nx_reference = (full ? 512 : 256) * scale;

    size_t npml = 24 * scale;
    size_t i_source = 36 * scale;
    size_t i_probe1 = 52 * scale;

    xt_sim reference = {0};
    xt_sim test = {0};

    int fail = 1;

    if (xt_create(
            &reference,
            nx_reference,
            ny,
            1,
            1e-3 / scale,
            6,
            npml) != 0) {
        goto cleanup;
    }

    if (xt_create(
            &test,
            nx_test,
            ny,
            1,
            1e-3 / scale,
            6,
            npml) != 0) {
        goto cleanup;
    }

    double courant =
        C0 * reference.g->dt / reference.g->dx;

    double ky =
        requested_angle == 0.0
        ? 0.0
        : 2.0 * XT_PI / ny;

    double target_kx =
        requested_angle == 0.0
        ? 2.0 * XT_PI / (32.0 * scale)
        : ky / tan(requested_angle * XT_PI / 180.0);

    double target_omega_dt = 2.0 * asin(
        courant * sqrt(
            pow(sin(target_kx / 2.0), 2.0)
            + pow(sin(ky / 2.0), 2.0)
        )
    );

    if (!isfinite(target_omega_dt) || target_omega_dt <= 0.0) {
        goto cleanup;
    }

    size_t period = (size_t)llround(
        2.0 * XT_PI / target_omega_dt
    );

    double omega_dt =
        requested_angle == 0.0
        ? 2.0 * XT_PI / period
        : target_omega_dt;

    double kx = xt_kx(omega_dt, courant, 1.0, ky);

    if (!isfinite(kx)) {
        fprintf(stderr, "Requested mode is not propagating on this grid\n");
        goto cleanup;
    }

    double actual_angle =
        atan2(ky, kx) * 180.0 / XT_PI;

    size_t probe_separation = (size_t)llround(0.7 / kx);

    if (probe_separation < 1) {
        probe_separation = 1;
    }

    size_t i_probe2 = i_probe1 + probe_separation;

    if (i_probe2 >= nx_test - npml - 3) {
        fprintf(stderr, "Insufficient space for oblique mode probes\n");
        goto cleanup;
    }

    double group_velocity_x =
        courant * courant * sin(kx) / sin(omega_dt);

    if (!isfinite(group_velocity_x) || group_velocity_x <= 0.0) {
        goto cleanup;
    }

    size_t ramp_steps = 8 * period;

    size_t settle_steps = ramp_steps + (size_t)ceil(
        6.0 * nx_reference / group_velocity_x
    );

    size_t window_steps = 12 * period;
    size_t total_steps = settle_steps + 2 * window_steps;

    double complex sums[2][2][2] = {{{0}}};

    fdtd_grid_t *grids[2] = {
        reference.g,
        test.g
    };

    size_t probe_positions[2] = {
        i_probe1,
        i_probe2
    };

    for (size_t n = 0; n < total_steps; n++) {
        xt_step(&reference);
        xt_step(&test);

        double phase = omega_dt * (n + 1);
        double envelope = 1.0;

        if (n < ramp_steps) {
            envelope = pow(
                sin(0.5 * XT_PI * (n + 1) / ramp_steps),
                2.0
            );
        }

        for (size_t j = 0; j < ny; j++) {
            double y_position = j + (tm ? 0.5 : 0.0);

            double source_value =
                envelope * cos(phase - ky * y_position);

            xt_e(reference.g, tm)[
                fdtd_index(reference.g, i_source, j, 0)
            ] += source_value;

            xt_e(test.g, tm)[
                fdtd_index(test.g, i_source, j, 0)
            ] += source_value;
        }

        if (n >= settle_steps) {
            size_t half = (n - settle_steps) / window_steps;

            double complex temporal_weight =
                cexp(I * phase) / window_steps;

            for (size_t g = 0; g < 2; g++) {
                for (size_t q = 0; q < 2; q++) {
                    double complex value = project_field(
                        grids[g],
                        probe_positions[q],
                        ky,
                        tm
                    );

                    if (!xt_finite_complex(value)) {
                        fprintf(stderr, "Nonfinite monitored field at step %zu\n", n);
                        goto cleanup;
                    }

                    sums[half][g][q] += value * temporal_weight;
                }
            }
        }
    }

    double complex forward[2][2];
    double complex backward[2][2];

    for (size_t h = 0; h < 2; h++) {
        for (size_t g = 0; g < 2; g++) {
            xt_modes(
                sums[h][g][0],
                sums[h][g][1],
                kx * probe_separation,
                &forward[h][g],
                &backward[h][g]
            );

            if (!xt_finite_complex(forward[h][g]) ||
                !xt_finite_complex(backward[h][g])) {
                goto cleanup;
            }
        }
    }

    if (cabs(forward[1][0]) < 1e-12 ||
        cabs(forward[1][1]) < 1e-12) {
        fprintf(stderr, "Incident mode amplitude is too small\n");
        goto cleanup;
    }

    double reference_reflection =
        cabs(backward[1][0] / forward[1][0]);

    double test_reflection =
        cabs(backward[1][1] / forward[1][1]);

    double backward_difference =
        cabs(
            (backward[1][1] - backward[1][0])
            / forward[1][0]
        );

    double window_drift = 0.0;

    for (size_t g = 0; g < 2; g++) {
        double drift = fmax(
            cabs(backward[1][g] - backward[0][g]),
            cabs(forward[1][g] - forward[0][g])
        ) / cabs(forward[1][g]);

        if (drift > window_drift) {
            window_drift = drift;
        }
    }

    fail = !(
        isfinite(test_reflection) &&
        isfinite(reference_reflection) &&
        isfinite(backward_difference) &&
        isfinite(window_drift) &&
        test_reflection <= 0.03 &&
        reference_reflection <= 0.03 &&
        window_drift <= 0.005 &&
        fabs(actual_angle - requested_angle) <= 1.5
    );

    printf(
        "%s requested %.1f actual %.4f deg: |r| %.5g (%.2f dB), reference %.5g, drift %.3g, steps %zu %s\n",
        tm ? "TM" : "TE",
        requested_angle,
        actual_angle,
        test_reflection,
        xt_db(test_reflection),
        reference_reflection,
        window_drift,
        total_steps,
        fail ? "FAIL" : "PASS"
    );

    if (window_drift > 0.1 * test_reflection) {
        printf("  Note: window drift is not small relative to reflection; do not treat this as a precision reflection measurement.\n");
    }

    if (fprintf(
            csv,
            "%s,%.17g,%.17g,%zu,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d\n",
            tm ? "TM" : "TE",
            requested_angle,
            actual_angle,
            scale,
            total_steps,
            omega_dt / (2.0 * XT_PI * reference.g->dt),
            test_reflection,
            xt_db(test_reflection),
            reference_reflection,
            backward_difference,
            window_drift,
            !fail) < 0) {
        fail = 1;
    }

cleanup:
    if (fail) {
        fprintf(stderr, "Oblique case failed or could not complete: %s %.1f degrees\n", tm ? "TM" : "TE", requested_angle);
    }

    xt_destroy(&reference);
    xt_destroy(&test);

    return fail;
}

int main(int argc, char **argv)
{
    int full = 0;

    if (xt_parse_full(argc, argv, &full) != 0) {
        fprintf(stderr, "Usage: %s [--full]\n", argv[0]);
        return 2;
    }

    FILE *csv = fopen("cpml_oblique_summary.csv", "w");

    if (!csv) {
        perror("cpml_oblique_summary.csv");
        return 2;
    }

    int fail = fprintf(
        csv,
        "polarization,requested_angle_deg,numerical_wavevector_angle_deg,refinement,steps,frequency_hz,reflection_amplitude,reflection_db,reference_reflection_amplitude,backward_difference_normalized,window_drift,pass\n"
    ) < 0;

    const double quick_angles[] = {
        30.0,
        60.0
    };

    const double full_angles[] = {
        0.0,
        15.0,
        30.0,
        45.0,
        60.0,
        75.0,
        85.0
    };

    size_t count = full ? 7 : 2;

    for (size_t q = 0; q < count; q++) {
        double angle = full ? full_angles[q] : quick_angles[q];

        for (int tm = 0; tm < 2; tm++) {
            fail |= run_case(csv, angle, tm, full);
        }
    }

    if (fclose(csv) != 0) {
        perror("Closing oblique summary");
        fail = 1;
    }

    return fail ? 1 : 0;
}