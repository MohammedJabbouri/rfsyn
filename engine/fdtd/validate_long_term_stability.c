#include "extended_test.h"

static int run_case(
    FILE * summary,
    const char * name,
    int closed,
    double eps_r,
    int has_interface,
    size_t total_steps,
    int full) {
	xt_sim s = {
		0
	};

	FILE * csv = NULL;
	double * epsilon = NULL;

	int fail = 1;

	size_t nx = closed ? 8 : (full ? 128 : 64);
	size_t ny = closed ? 8 : 1;
	size_t nz = closed ? 8 : 1;

	size_t npml = closed ? 0 : 16;
	unsigned periodic_axes = closed ? 7u : 6u;

	if (xt_create( &
	               s,
	               nx,
	               ny,
	               nz,
	               1e-3,
	               periodic_axes,
	               npml) != 0) {
		goto cleanup;
	}

	if (has_interface) {
		if (xt_interface( & s, nx / 2, 1.0, eps_r) != 0) {
			goto cleanup;
		}
	} else {
		if (xt_medium( & s, eps_r) != 0) {
			goto cleanup;
		}
	}

	epsilon = malloc(nx * ny * nz * 3 * sizeof( * epsilon));

	if (!epsilon) {
		fprintf(stderr, "Energy coefficient allocation failed\n");
		goto cleanup;
	}

	for (size_t i = 0; i < nx; i++) {
		for (size_t j = 0; j < ny; j++) {
			for (size_t k = 0; k < nz; k++) {
				for (int c = 0; c < 3; c++) {
					material_coeff_t coefficient;

					if (materials_get_coefficients(
					            s.m,
					            s.g,
					            (field_component_t) c,
					            i,
					            j,
					            k, &
					            coefficient) != 0) {
						goto cleanup;
					}

					size_t z = (((i * ny + j) * nz + k) * 3) + c;

					epsilon[z] =
					    s.g -> dt / (s.g -> dx * coefficient.Cb);
				}
			}
		}
	}

	char filename[128];

	snprintf(filename, sizeof(filename), "stability_%s.csv", name);

	csv = fopen(filename, "w");

	if (!csv) {
		perror(filename);
		goto cleanup;
	}

	if (fprintf(
	            csv,
	            "step,t,max_abs_e,max_abs_h,diagnostic_energy,interior_energy,source_active\n"
	        ) < 0) {
		goto cleanup;
	}

	double * electric[3] = {
		s.g -> ex,
		s.g -> ey,
		s.g -> ez
	};

	double * magnetic[3] = {
		s.g -> hx,
		s.g -> hy,
		s.g -> hz
	};

	uint32_t seed = 12345u;

	if (closed) {
		for (size_t i = 0; i < nx; i++) {
			for (size_t j = 0; j < ny; j++) {
				for (size_t k = 0; k < nz; k++) {
					size_t idx = fdtd_index(s.g, i, j, k);

					for (size_t c = 0; c < 3; c++) {
						seed = 1664525u * seed + 1013904223u;

						electric[c][idx] =
						    ((double)(seed >> 8) / 16777216.0 - 0.5) *
						    0.1;
					}
				}
			}
		}
	}

	double maximum_energy = 0.0;
	double reference_energy = 0.0;
	double late_peak_energy = 0.0;
	double final_energy = 0.0;

	double previous_window_peak = 0.0;
	double current_window_peak = 0.0;

	size_t consecutive_growth = 0;
	size_t maximum_growth_run = 0;
	size_t completed_windows = 0;

	size_t sample_stride = closed ? 50 : 100;
	size_t source_end = closed ? 0 : 400;

	for (size_t n = 0; n < total_steps; n++) {
		xt_step( & s);

		if (!closed && n < source_end) {
			double u = ((double) n - 150.0) / 30.0;
			double source_value = -2.0 * u * exp(-u * u);

			size_t idx = fdtd_index(s.g, nx / 2 - 4, 0, 0);
			s.g -> ez[idx] += source_value;
		}

		if (n % sample_stride != 0 && n + 1 != total_steps) {
			continue;
		}

		double energy = 0.0;
		double interior_energy = 0.0;
		double maximum_e = 0.0;
		double maximum_h = 0.0;

		for (size_t i = 0; i < nx; i++) {
			for (size_t j = 0; j < ny; j++) {
				for (size_t k = 0; k < nz; k++) {
					size_t idx = fdtd_index(s.g, i, j, k);
					size_t z = ((i * ny + j) * nz + k) * 3;

					double cell_energy = 0.0;

					for (size_t c = 0; c < 3; c++) {
						double ev = electric[c][idx];
						double hv = magnetic[c][idx];

						if (!isfinite(ev) || !isfinite(hv)) {
							fprintf(stderr, "Nonfinite field in %s at step %zu\n", name, n);
							goto cleanup;
						}

						if (fabs(ev) > maximum_e) {
							maximum_e = fabs(ev);
						}

						if (fabs(hv) > maximum_h) {
							maximum_h = fabs(hv);
						}

						cell_energy += 0.5 * (
						                   epsilon[z + c] * ev * ev +
						                   MU0 * hv * hv
						               );
					}

					energy += cell_energy;

					if (closed || (i >= npml && i < nx - npml)) {
						interior_energy += cell_energy;
					}
				}
			}
		}

		double cell_volume = pow(s.g -> dx, 3.0);

		energy *= cell_volume;
		interior_energy *= cell_volume;

		if (!isfinite(energy) || !isfinite(interior_energy)) {
			fprintf(stderr, "Nonfinite energy in %s at step %zu\n", name, n);
			goto cleanup;
		}

		if (energy > maximum_energy) {
			maximum_energy = energy;
		}

		if (closed && n == 0) {
			reference_energy = energy;
		}

		if (!closed && n <= source_end + 1000) {
			reference_energy = maximum_energy;
		}

		if (n >= total_steps / 2 && energy > late_peak_energy) {
			late_peak_energy = energy;
		}

		if (energy > current_window_peak) {
			current_window_peak = energy;
		}

		if (n > 0 && n % 1000 == 0) {
			/*
			 * Require sustained growth above a diagnostic floor.
			 * Below-floor beating alone is not classified as instability.
			 */
			if (n >= total_steps / 2 &&
			        previous_window_peak > reference_energy * 1e-8 &&
			        current_window_peak > previous_window_peak * 1.2) {
				consecutive_growth++;
			} else {
				consecutive_growth = 0;
			}

			if (consecutive_growth > maximum_growth_run) {
				maximum_growth_run = consecutive_growth;
			}

			previous_window_peak = current_window_peak;
			current_window_peak = 0.0;
			completed_windows++;
		}

		final_energy = energy;

		if (fprintf(
		            csv,
		            "%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%d\n",
		            n,
		            (n + 1) * s.g -> dt,
		            maximum_e,
		            maximum_h,
		            energy,
		            interior_energy,
		            !closed && n < source_end) < 0) {
			goto cleanup;
		}

		if (reference_energy > 0.0 &&
		        energy > reference_energy * 100.0) {
			fprintf(stderr, "Large energy growth detected in %s at step %zu\n", name, n);
			goto cleanup;
		}
	}

	if (!(reference_energy > 0.0) || completed_windows == 0) {
		fprintf(stderr, "Insufficient nonzero reference data in %s\n", name);
		goto cleanup;
	}

	double maximum_ratio = maximum_energy / reference_energy;
	double late_ratio = late_peak_energy / reference_energy;
	double final_ratio = final_energy / reference_energy;

	if (closed) {
		fail = maximum_ratio > 1.5 || final_ratio < 0.5;
	} else {
		fail =
		    late_ratio > 1e-4 ||
		    final_ratio > 1e-6 ||
		    maximum_growth_run >= 5;
	}

	printf(
	    "%s steps=%zu max/reference=%.6g late_peak/reference=%.6g final/reference=%.6g growth_windows=%zu %s\n",
	    name,
	    total_steps,
	    maximum_ratio,
	    late_ratio,
	    final_ratio,
	    maximum_growth_run,
	    fail ? "FAIL" : "PASS"
	);

	if (fprintf(
	            summary,
	            "%s,%zu,%.17g,%.17g,%.17g,%zu,%d\n",
	            name,
	            total_steps,
	            maximum_ratio,
	            late_ratio,
	            final_ratio,
	            maximum_growth_run,
	            !fail) < 0) {
		fail = 1;
	}

cleanup:
	if (csv && fclose(csv) != 0) {
		perror("Closing stability history");
		fail = 1;
	}

	free(epsilon);
	xt_destroy( & s);

	if (fail) {
		fprintf(stderr, "Stability case failed or could not complete: %s\n", name);
	}

	return fail;
}

int main(int argc, char ** argv) {
	size_t steps = 100000;
	int full = 0;

	if (xt_parse_stability(argc, argv, & steps, & full) != 0) {
		fprintf(stderr, "Usage: %s [--steps 100000] [--full]\n", argv[0]);
		return 2;
	}

	FILE * summary = fopen("long_term_stability_summary.csv", "w");

	if (!summary) {
		perror("long_term_stability_summary.csv");
		return 2;
	}

	int fail = fprintf(
	               summary,
	               "case,steps,max_energy_over_reference,late_peak_over_reference,final_over_reference,growth_windows,pass\n"
	           ) < 0;

	fail |= run_case(
	            summary,
	            "open_vacuum",
	            0,
	            1.0,
	            0,
	            steps,
	            full
	        );

	fail |= run_case(
	            summary,
	            "open_eps25",
	            0,
	            25.0,
	            0,
	            steps,
	            full
	        );

	fail |= run_case(
	            summary,
	            "open_interface25",
	            0,
	            25.0,
	            1,
	            steps,
	            full
	        );

	fail |= run_case(
	            summary,
	            "periodic_vacuum_3d",
	            1,
	            1.0,
	            0,
	            steps,
	            full
	        );

	fail |= run_case(
	            summary,
	            "periodic_interface25_3d",
	            1,
	            25.0,
	            1,
	            steps,
	            full
	        );

	if (fclose(summary) != 0) {
		perror("Closing stability summary");
		fail = 1;
	}

	return fail ? 1 : 0;
}