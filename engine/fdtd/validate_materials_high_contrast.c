#include "extended_test.h"

static int run_case(FILE *csv, double eps_left, double eps_right)
{
    const size_t scale = RFSYN_VALIDATION_REFINEMENT;

    const size_t nx = 256 * scale;
    const size_t npml = 24 * scale;

    const size_t i_source = 40 * scale;
    const size_t i_interface = 112 * scale;

    const size_t i_probe1 = 64 * scale;
    const size_t i_probe2 = 80 * scale;

    const size_t i_trans1 = 144 * scale;
    const size_t i_trans2 = 160 * scale;

    xt_sim s = {0};
    int fail = 1;

    if (xt_create(&s, nx, 1, 1, 1e-3 / scale, 6, npml) != 0) {
        goto cleanup;
    }

    if (xt_interface(&s, i_interface, eps_left, eps_right) != 0) {
        goto cleanup;
    }

    double courant = C0 * s.g->dt / s.g->dx;
    double maximum_index = sqrt(fmax(eps_left, eps_right));

    size_t period = (size_t)ceil(
        40.0 * maximum_index * scale / courant
    );

    double omega_dt = 2.0 * XT_PI / period;

    double k_left = xt_kx(omega_dt, courant, eps_left, 0.0);
    double k_right = xt_kx(omega_dt, courant, eps_right, 0.0);

    if (!isfinite(k_left) || !isfinite(k_right)) {
        goto cleanup;
    }

    if (fabs(sin(k_left * (i_probe2 - i_probe1))) < 0.02 ||
        fabs(sin(k_right * (i_trans2 - i_trans1))) < 0.02) {
        fprintf(stderr, "Probe spacing is ill-conditioned for this frequency\n");
        goto cleanup;
    }

    size_t ramp_steps = 8 * period;

    size_t settle_steps = ramp_steps + (size_t)ceil(
        8.0 * nx * maximum_index / courant
    );

    size_t window_steps = 8 * period;
    size_t total_steps = settle_steps + 2 * window_steps;

    double complex sums[2][4] = {{0}};

    const size_t probe_positions[4] = {
        i_probe1,
        i_probe2,
        i_trans1,
        i_trans2
    };

    for (size_t n = 0; n < total_steps; n++) {
        xt_step(&s);

        double phase = omega_dt * (n + 1);
        double envelope = 1.0;

        if (n < ramp_steps) {
            envelope = pow(
                sin(0.5 * XT_PI * (n + 1) / ramp_steps),
                2.0
            );
        }

        size_t source_index = fdtd_index(s.g, i_source, 0, 0);
        s.g->ez[source_index] += envelope * cos(phase);

        if (n >= settle_steps) {
            size_t half = (n - settle_steps) / window_steps;
            double complex temporal_weight =
                cexp(I * phase) / window_steps;

            for (size_t q = 0; q < 4; q++) {
                size_t idx = fdtd_index(
                    s.g,
                    probe_positions[q],
                    0,
                    0
                );

                double value = s.g->ez[idx];

                if (!isfinite(value)) {
                    fprintf(stderr, "Nonfinite monitored field at step %zu\n", n);
                    goto cleanup;
                }

                sums[half][q] += value * temporal_weight;
            }
        }
    }

    double complex forward_left[2];
    double complex backward_left[2];
    double complex forward_right[2];
    double complex backward_right[2];

    for (size_t h = 0; h < 2; h++) {
        xt_modes(
            sums[h][0],
            sums[h][1],
            k_left * (i_probe2 - i_probe1),
            &forward_left[h],
            &backward_left[h]
        );

        xt_modes(
            sums[h][2],
            sums[h][3],
            k_right * (i_trans2 - i_trans1),
            &forward_right[h],
            &backward_right[h]
        );

        if (!xt_finite_complex(forward_left[h]) ||
            !xt_finite_complex(backward_left[h]) ||
            !xt_finite_complex(forward_right[h]) ||
            !xt_finite_complex(backward_right[h])) {
            goto cleanup;
        }
    }

    if (cabs(forward_left[1]) < 1e-12) {
        fprintf(stderr, "Incident amplitude is too small\n");
        goto cleanup;
    }

    double complex reflection =
        backward_left[1] / forward_left[1]
        * cexp(-2.0 * I * k_left * (i_interface - i_probe1));

    double transmission_amplitude =
        cabs(forward_right[1] / forward_left[1]);

    double R = pow(cabs(reflection), 2.0);

    double T =
        sqrt(eps_right / eps_left)
        * transmission_amplitude
        * transmission_amplitude;

    double r_theory =
        (sqrt(eps_left) - sqrt(eps_right))
        / (sqrt(eps_left) + sqrt(eps_right));

    double R_theory = r_theory * r_theory;
    double T_theory = 1.0 - R_theory;

    double amplitude_error =
        fabs(cabs(reflection) - fabs(r_theory));

    double power_error = fmax(
        fabs(R - R_theory),
        fabs(T - T_theory)
    );

    double window_drift = fmax(
        cabs(forward_left[1] - forward_left[0]),
        cabs(backward_left[1] - backward_left[0])
    ) / cabs(forward_left[1]);

    double transmitted_backward_ratio =
        cabs(backward_right[1])
        / fmax(cabs(forward_right[1]), 1e-30);

    double phase_error_deg = NAN;

    if (r_theory != 0.0) {
        phase_error_deg =
            fabs(carg(reflection / r_theory))
            * 180.0 / XT_PI;
    }

    fail = !(
        isfinite(R) &&
        isfinite(T) &&
        isfinite(window_drift) &&
        isfinite(transmitted_backward_ratio) &&
        amplitude_error <= 0.01 &&
        power_error <= 0.02 &&
        fabs(R + T - 1.0) <= 0.02 &&
        window_drift <= 0.002 &&
        transmitted_backward_ratio <= 0.01 &&
        (r_theory == 0.0 || phase_error_deg <= 3.0)
    );

    printf(
        "eps %.0f -> %.0f: |r| %.8g theory %.8g R %.8g T %.8g drift %.3g phase_error %.4g deg %s\n",
        eps_left,
        eps_right,
        cabs(reflection),
        fabs(r_theory),
        R,
        T,
        window_drift,
        phase_error_deg,
        fail ? "FAIL" : "PASS"
    );

    if (fprintf(csv, "%.17g,%.17g,%zu,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d\n", eps_left, eps_right, scale, total_steps, 1.0 / (period * s.g->dt), cabs(reflection), fabs(r_theory), R, T, R_theory, phase_error_deg, window_drift, transmitted_backward_ratio, !fail) < 0) {fail = 1;
    }

cleanup:
    if (fail) {
        fprintf(stderr, "High-contrast case failed or could not complete: eps %.0f -> %.0f\n", eps_left, eps_right);
    }

    xt_destroy(&s);
    return fail;
}

int main(int argc, char **argv)
{
    int full = 0;

    if (xt_parse_full(argc, argv, &full) != 0) {
        fprintf(stderr, "Usage: %s [--full]\n", argv[0]);
        return 2;
    }

    FILE *csv = fopen("materials_high_contrast_summary.csv", "w");

    if (!csv) {
        perror("materials_high_contrast_summary.csv");
        return 2;
    }

    int fail = fprintf(csv, "eps_left,eps_right,refinement,steps,frequency_hz,r_magnitude,r_theory_magnitude,R,T,R_theory,reflection_phase_error_deg,window_drift,backward_transmitted_ratio,pass\n" ) < 0;

    const double contrasts[] = {
        1.0,
        4.0,
        10.0,
        25.0,
        50.0,
        100.0
    };

    size_t count = full ? 6 : 4;

    for (size_t q = 0; q < count; q++) {
        fail |= run_case(csv, 1.0, contrasts[q]);

        if (contrasts[q] != 1.0) {
            fail |= run_case(csv, contrasts[q], 1.0);
        }
    }

    if (fclose(csv) != 0) {
        perror("Closing high-contrast summary");
        fail = 1;
    }

    return fail ? 1 : 0;
}