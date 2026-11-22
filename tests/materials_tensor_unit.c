#include "../engine/fdtd/extended_test.h"

#include <float.h>

#include <complex.h>
#include "../engine/core/constants.h"

static size_t checks = 0;
static size_t failures = 0;

#define CHECK(expr)                                                   \
    do {                                                              \
        checks++;                                                     \
        if (!(expr)) {                                                \
            fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #expr);    \
            failures++;                                               \
        }                                                             \
    } while (0)

#define REQUIRE(expr)                                                 \
    do {                                                              \
        checks++;                                                     \
        if (!(expr)) {                                                \
            fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #expr);    \
            failures++;                                               \
            goto cleanup;                                             \
        }                                                             \
    } while (0)

static void check_close(
    double actual,
    double expected,
    double atol,
    double rtol,
    int line) {
    checks++;

    if (!isfinite(actual) ||
        !isfinite(expected) ||
        fabs(actual - expected) >
        atol + rtol * fabs(expected)) {
        fprintf(
            stderr,
            "FAIL line %d: actual=%.17g expected=%.17g "
            "atol=%.3g rtol=%.3g\n",
            line, actual, expected, atol, rtol);
        failures++;
    }
}

#define CLOSE(actual, expected, atol, rtol) \
    check_close((actual), (expected), (atol), (rtol), __LINE__)

static size_t grid_storage(const fdtd_grid_t *g) {
    return (g->nx + 1) * (g->ny + 1) * (g->nz + 1);
}

/*
 * All fixtures use periodic z, with one active z plane.
 * Ex/Ey at k=0 therefore really are updated.
 *
 * The test edges lie well inside any CPML slabs.
 */
static int create_fixture(xt_sim *s, int with_cpml) {
    return xt_create(
        s,
        24, 24, 1,
        0.001,
        4u,
        with_cpml ? 4u : 0u);
}

static void update_e(xt_sim *s) {
    if (s->p) {
        fdtd_update_e_cpml_materials(s->g, s->p, s->m);
    } else {
        fdtd_update_e_materials(s->g, s->m);
    }
}

static void compare_fields(
    const fdtd_grid_t *a,
    const fdtd_grid_t *b,
    double atol,
    double rtol) {
    const double *af[6] = {
        a->ex, a->ey, a->ez, a->hx, a->hy, a->hz
    };

    const double *bf[6] = {
        b->ex, b->ey, b->ez, b->hx, b->hy, b->hz
    };

    static const char *names[6] = {
        "Ex", "Ey", "Ez", "Hx", "Hy", "Hz"
    };

    size_t total = grid_storage(a);

    for (size_t f = 0; f < 6; f++) {
        int valid = 1;
        size_t bad_index = 0;
        double max_error = 0.0;

        for (size_t q = 0; q < total; q++) {
            double av = af[f][q];
            double bv = bf[f][q];

            if (!isfinite(av) || !isfinite(bv)) {
                valid = 0;
                bad_index = q;
                break;
            }

            double error = fabs(av - bv);

            if (error > max_error) {
                max_error = error;
            }

            double magnitude = fmax(fabs(av), fabs(bv));

            if (error > atol + rtol * magnitude) {
                valid = 0;
                bad_index = q;
                break;
            }
        }

        checks++;

        if (!valid) {
            fprintf(
                stderr,
                "FAIL field %s index %zu: %.17g vs %.17g "
                "(largest error examined %.17g)\n",
                names[f],
                bad_index,
                af[f][bad_index],
                bf[f][bad_index],
                max_error);
            failures++;
        }
    }
}

static void test_api(void) {
    xt_sim s = {0};
    material_tensor_t out = {0};

    REQUIRE(create_fixture(&s, 0) == 0);

    CHECK(materials_get_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0, &out) == 1);

    REQUIRE(materials_set_component_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0,
        0.6, 0.05, 0.5) == 0);

    REQUIRE(materials_get_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0, &out) == 0);

    CLOSE(out.m11, 0.6, 0.0, 0.0);
    CLOSE(out.m12, 0.05, 0.0, 0.0);
    CLOSE(out.m22, 0.5, 0.0, 0.0);

    CLOSE(
        out.scale,
        s.g->dt / (EPS0 * s.g->dx),
        0.0,
        16.0 * DBL_EPSILON);

    CHECK(materials_set_component_tensor(
        s.m, s.g, FIELD_EZ, 10, 10, 0,
        0.6, 0.0, 0.5) == -1);

    CHECK(materials_set_component_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0,
        0.3, 0.4, 0.3) == -1);

    CHECK(materials_set_component_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0,
        0.0, 0.0, 0.3) == -1);

    CHECK(materials_set_component_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0,
        -0.3, 0.0, 0.3) == -1);

    CHECK(materials_set_component_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0,
        0.3, NAN, 0.3) == -1);

    CHECK(materials_set_component_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0,
        INFINITY, 0.0, 0.3) == -1);

    CHECK(materials_set_component_tensor(
        s.m, s.g, FIELD_EX, 0, 10, 0,
        0.6, 0.0, 0.5) == -1);

    CHECK(materials_set_component_tensor(
        s.m, s.g, FIELD_EY, 10, 0, 0,
        0.6, 0.0, 0.5) == -1);

    CHECK(materials_set_component_tensor(
        s.m, s.g, FIELD_EX, s.g->nx, 10, 0,
        0.6, 0.0, 0.5) == -1);

    CHECK(materials_get_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0, NULL) == -1);

    /* Rejected writes must not change the existing entry. */
    REQUIRE(materials_get_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0, &out) == 0);

    CLOSE(out.m11, 0.6, 0.0, 0.0);
    CLOSE(out.m12, 0.05, 0.0, 0.0);
    CLOSE(out.m22, 0.5, 0.0, 0.0);

    /* Insert out of order, including two components at one index. */
    REQUIRE(materials_set_component_tensor(
        s.m, s.g, FIELD_EY, 8, 8, 0,
        0.7, -0.04, 0.4) == 0);

    REQUIRE(materials_set_component_tensor(
        s.m, s.g, FIELD_EX, 12, 12, 0,
        0.8, 0.02, 0.6) == 0);

    REQUIRE(materials_set_component_tensor(
        s.m, s.g, FIELD_EY, 10, 10, 0,
        0.9, 0.01, 0.4) == 0);

    REQUIRE(materials_set_component_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0,
        0.4, 0.0, 0.4) == 0);

    REQUIRE(materials_get_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0, &out) == 0);

    CLOSE(out.m11, 0.4, 0.0, 0.0);

    REQUIRE(materials_get_tensor(
        s.m, s.g, FIELD_EY, 10, 10, 0, &out) == 0);

    CLOSE(out.m11, 0.9, 0.0, 0.0);

    CHECK(materials_clear_component_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0) == 0);

    CHECK(materials_get_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0, &out) == 1);

    /* Clearing one component must preserve the other. */
    CHECK(materials_get_tensor(
        s.m, s.g, FIELD_EY, 10, 10, 0, &out) == 0);

    CHECK(materials_clear_component_tensor(
        s.m, s.g, FIELD_EX, 10, 10, 0) == 0);

    REQUIRE(materials_set_box(
        s.m, s.g,
        9, 11, 9, 11, 0, 1,
        4.0, 0.0) == 0);

    CHECK(materials_get_tensor(
        s.m, s.g, FIELD_EY, 10, 10, 0, &out) == 1);

    REQUIRE(materials_get_tensor(
        s.m, s.g, FIELD_EY, 8, 8, 0, &out) == 0);

    CLOSE(out.m11, 0.7, 0.0, 0.0);

    REQUIRE(materials_get_tensor(
        s.m, s.g, FIELD_EX, 12, 12, 0, &out) == 0);

    CLOSE(out.m11, 0.8, 0.0, 0.0);

cleanup:
    xt_destroy(&s);
}

/*
 * Hz(i,j) = A*i + B*j + Q*i^2 + R*j^2.
 *
 * Analytic centered derivatives:
 *   dHz/di = A + 2*Q*i
 *   dHz/dj = B + 2*R*j
 *
 * Yee backward differences:
 *   Hz(i,j) - Hz(i,j-1) = B + R*(2*j-1)
 *   Hz(i,j) - Hz(i-1,j) = A + Q*(2*i-1)
 */
static void test_analytic_update(int with_cpml, int quadratic) {
    xt_sim s = {0};

    const size_t i = 10;
    const size_t j = 11;

    const double A = 0.7;
    const double B = -0.2;
    const double Q = quadratic ? 0.03 : 0.0;
    const double R = quadratic ? -0.02 : 0.0;

    const double mxx = 0.65;
    const double mxy = 0.12;
    const double myy = 0.40;

    const double ex_old = 0.13;
    const double ey_old = -0.17;

    REQUIRE(create_fixture(&s, with_cpml) == 0);

    for (size_t ii = 0; ii <= s.g->nx; ii++) {
        for (size_t jj = 0; jj <= s.g->ny; jj++) {
            for (size_t kk = 0; kk <= s.g->nz; kk++) {
                double x = (double) ii;
                double y = (double) jj;

                s.g->hz[fdtd_index(s.g, ii, jj, kk)] =
                        A * x + B * y + Q * x * x + R * y * y;
            }
        }
    }

    size_t idx = fdtd_index(s.g, i, j, 0);

    s.g->ex[idx] = ex_old;
    s.g->ey[idx] = ey_old;

    /* Uniform patch supplies all four reciprocal neighbors. */
    for (size_t ii = i - 1; ii <= i + 1; ii++) {
        for (size_t jj = j - 1; jj <= j + 1; jj++) {
            REQUIRE(materials_set_component_tensor(
                s.m, s.g, FIELD_EX, ii, jj, 0,
                mxx, mxy, myy) == 0);
            REQUIRE(materials_set_component_tensor(
                s.m, s.g, FIELD_EY, ii, jj, 0,
                myy, mxy, mxx) == 0);
        }
    }

    double scale = s.g->dt / (EPS0 * s.g->dx);

    double ex_aligned = B + R * (2.0 * j - 1.0);
    double ex_cross = -(A + 2.0 * Q * i);

    double ey_aligned = -(A + Q * (2.0 * i - 1.0));
    double ey_cross = B + 2.0 * R * j;

    double expected_ex =
            ex_old + scale * (mxx * ex_aligned + mxy * ex_cross);

    double expected_ey =
            ey_old + scale * (myy * ey_aligned + mxy * ey_cross);

    update_e(&s);

    fprintf(
        stdout,
        "analytic update: cpml=%d quadratic=%d\n",
        with_cpml, quadratic);

    CLOSE(s.g->ex[idx], expected_ex, 1e-10, 1e-12);
    CLOSE(s.g->ey[idx], expected_ey, 1e-10, 1e-12);

    /* Neighboring scalar edges must still take the scalar path. */
    size_t ni = i + 2;
    size_t nj = j + 2;
    size_t neighbor = fdtd_index(s.g, ni, nj, 0);

    CLOSE(
        s.g->ex[neighbor],
        scale * (B + R * (2.0 * nj - 1.0)),
        1e-10,
        1e-12);

    CLOSE(
        s.g->ey[neighbor],
        -scale * (A + Q * (2.0 * ni - 1.0)),
        1e-10,
        1e-12);

cleanup:
    xt_destroy(&s);
}

/*
 * Compare two otherwise identical simulations:
 *   a: scalar epsilon everywhere;
 *   b: same scalar epsilon, with equivalent diagonal tensor overrides.
 *
 * This tests diagonal reduction, not the accuracy of interface mixing.
 */
static void test_diagonal_invariance(int with_cpml) {
    xt_sim a = {0};
    xt_sim b = {0};

    const double eps_r = 2.5;
    const double inverse_eps = 1.0 / eps_r;

    REQUIRE(create_fixture(&a, with_cpml) == 0);
    REQUIRE(create_fixture(&b, with_cpml) == 0);

    REQUIRE(xt_medium(&a, eps_r) == 0);
    REQUIRE(xt_medium(&b, eps_r) == 0);

    for (size_t i = 8; i < 16; i++) {
        for (size_t j = 8; j < 16; j++) {
            REQUIRE(materials_set_component_tensor(
                b.m, b.g, FIELD_EX, i, j, 0,
                inverse_eps, 0.0, inverse_eps) == 0);

            REQUIRE(materials_set_component_tensor(
                b.m, b.g, FIELD_EY, i, j, 0,
                inverse_eps, 0.0, inverse_eps) == 0);
        }
    }

    /* Smooth, localized initial Hz field, identical in both grids. */
    for (size_t i = 0; i <= a.g->nx; i++) {
        for (size_t j = 0; j <= a.g->ny; j++) {
            double x = ((double) i - 12.0) / 3.0;
            double y = ((double) j - 12.0) / 3.0;
            double value = exp(-(x * x + y * y));

            for (size_t k = 0; k <= a.g->nz; k++) {
                size_t idx = fdtd_index(a.g, i, j, k);
                a.g->hz[idx] = value;
                b.g->hz[idx] = value;
            }
        }
    }

    for (int step = 0; step < 80; step++) {
        xt_step(&a);
        xt_step(&b);

        compare_fields(a.g, b.g, 1e-10, 1e-11);
    }

    fprintf(
        stdout,
        "diagonal invariance: cpml=%d steps=80\n",
        with_cpml);

cleanup:
    xt_destroy(&a);
    xt_destroy(&b);
}

/*
 * Periodic-seam oracle.
 *
 * Hz(i,j) =
 *     sin(ax*i + ay*j + phase_a)
 *   + amplitude_b*cos(bx*i + by*j + phase_b)
 *
 * Frequencies contain integer numbers of periods over each extent.
 * Expected Yee differences are evaluated with trigonometric identities,
 * not with the production neighbor-index helpers.
 */
static double seam_hz_value(
    double i,
    double j,
    double ax,
    double ay,
    double bx,
    double by) {
    const double phase_a = 0.37;
    const double phase_b = -0.29;
    const double amplitude_b = 0.31;

    return
            sin(ax * i + ay * j + phase_a) +
            amplitude_b * cos(bx * i + by * j + phase_b);
}

static void seam_expected_curls(
    size_t i,
    size_t j,
    double ax,
    double ay,
    double bx,
    double by,
    double *ex_aligned,
    double *ex_cross,
    double *ey_aligned,
    double *ey_cross) {
    const double phase_a = 0.37;
    const double phase_b = -0.29;
    const double amplitude_b = 0.31;

    double x = (double) i;
    double y = (double) j;

    /*
     * Ex is at (i + 1/2, j).
     * The Hz samples in the stored-index convention give phases
     * centered at (i, j - 1/2) for these differences.
     */
    double ex_phase_a = ax * x + ay * (y - 0.5) + phase_a;
    double ex_phase_b = bx * x + by * (y - 0.5) + phase_b;

    *ex_aligned =
            2.0 * sin(0.5 * ay) * cos(ex_phase_a) -
            2.0 * amplitude_b * sin(0.5 * by) * sin(ex_phase_b);

    *ex_cross = -(
        sin(ax) * cos(0.5 * ay) * cos(ex_phase_a) -
        amplitude_b * sin(bx) * cos(0.5 * by) * sin(ex_phase_b)
    );

    /*
     * Ey aligned curl is the negative x-directed backward difference.
     * Its cross curl averages y-directed differences onto Ey.
     */
    double ey_phase_a = ax * (x - 0.5) + ay * y + phase_a;
    double ey_phase_b = bx * (x - 0.5) + by * y + phase_b;

    *ey_aligned = -(
        2.0 * sin(0.5 * ax) * cos(ey_phase_a) -
        2.0 * amplitude_b * sin(0.5 * bx) * sin(ey_phase_b)
    );

    *ey_cross =
            sin(ay) * cos(0.5 * ax) * cos(ey_phase_a) -
            amplitude_b * sin(by) * cos(0.5 * bx) * sin(ey_phase_b);
}

static void test_periodic_tensor_seams(
    unsigned periodic_axes,
    int with_cpml,
    double mxy) {
    xt_sim s = {0};

    const size_t nx = 24;
    const size_t ny = 24;

    const double mxx = 0.65;
    const double myy = 0.40;

    const double ex_old = 0.13;
    const double ey_old = -0.17;

    /*
     * Distinct, nonzero spatial frequencies make x/y swaps and
     * one-cell shifts observable.
     */
    const double ax = 2.0 * XT_PI * 2.0 / (double) nx;
    const double ay = 2.0 * XT_PI * 3.0 / (double) ny;
    const double bx = 2.0 * XT_PI * 5.0 / (double) nx;
    const double by = 2.0 * XT_PI * 4.0 / (double) ny;

    fprintf(
        stdout,
        "periodic tensor seams: axes=%u cpml=%d mxy=%+.3f\n",
        periodic_axes, with_cpml, mxy);
    fflush(stdout);

    REQUIRE(xt_create(
        &s,
        nx, ny, 1,
        0.001,
        periodic_axes,
        with_cpml ? 4u : 0u) == 0);

    /*
     * Canonical Hz samples are filled below. Leave all other storage
     * as NaN so accidental duplicate-endpoint/padding reads propagate
     * to the checked electric fields.
     */
    for (size_t q = 0; q < grid_storage(s.g); q++) {
        s.g->hz[q] = NAN;
    }

    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            s.g->hz[fdtd_index(s.g, i, j, 0)] =
                    seam_hz_value(
                        (double) i,
                        (double) j,
                        ax, ay, bx, by);
        }
    }

    /*
     * On a periodic axis, test every canonical coordinate.
     *
     * On a nonperiodic axis, use [6,18), safely away from the
     * four-cell CPML slabs and hard boundaries.
     */
    size_t ibegin = (periodic_axes & 1u) ? 0u : 6u;
    size_t iend = (periodic_axes & 1u) ? nx : 18u;
    size_t jbegin = (periodic_axes & 2u) ? 0u : 6u;
    size_t jend = (periodic_axes & 2u) ? ny : 18u;

    size_t ti0 = (periodic_axes & 1u) ? ibegin : ibegin - 1;
    size_t ti1 = (periodic_axes & 1u) ? iend : iend + 1;
    size_t tj0 = (periodic_axes & 2u) ? jbegin : jbegin - 1;
    size_t tj1 = (periodic_axes & 2u) ? jend : jend + 1;

    for (size_t i = ti0; i < ti1; i++) {
        for (size_t j = tj0; j < tj1; j++) {
            size_t idx = fdtd_index(s.g, i, j, 0);

            s.g->ex[idx] = ex_old;
            s.g->ey[idx] = ey_old;

            REQUIRE(materials_set_component_tensor(
                s.m, s.g, FIELD_EX, i, j, 0,
                mxx, mxy, myy) == 0);

            REQUIRE(materials_set_component_tensor(
                s.m, s.g, FIELD_EY, i, j, 0,
                myy, mxy, mxx) == 0);
        }
    }

    /*
     * Periodic duplicate endpoints must be rejected, not accepted
     * as additional tensor coordinates or silently aliased onto zero.
     */
    if (periodic_axes & 1u) {
        CHECK(materials_set_component_tensor(
            s.m, s.g, FIELD_EX, nx, jbegin, 0,
            0.9, 0.01, 0.8) == -1);

        CHECK(materials_set_component_tensor(
            s.m, s.g, FIELD_EY, nx, jbegin, 0,
            0.8, 0.01, 0.9) == -1);
    }

    if (periodic_axes & 2u) {
        CHECK(materials_set_component_tensor(
            s.m, s.g, FIELD_EX, ibegin, ny, 0,
            0.9, 0.01, 0.8) == -1);

        CHECK(materials_set_component_tensor(
            s.m, s.g, FIELD_EY, ibegin, ny, 0,
            0.8, 0.01, 0.9) == -1);
    }

    /*
     * Check that rejected duplicate-endpoint writes did not alter
     * the canonical entry.
     */
    {
        material_tensor_t out = {0};

        REQUIRE(materials_get_tensor(
            s.m, s.g, FIELD_EX,
            ibegin, jbegin, 0, &out) == 0);

        CLOSE(out.m11, mxx, 0.0, 0.0);
        CLOSE(out.m12, mxy, 0.0, 0.0);
        CLOSE(out.m22, myy, 0.0, 0.0);

        REQUIRE(materials_get_tensor(
            s.m, s.g, FIELD_EY,
            ibegin, jbegin, 0, &out) == 0);

        CLOSE(out.m11, myy, 0.0, 0.0);
        CLOSE(out.m12, mxy, 0.0, 0.0);
        CLOSE(out.m22, mxx, 0.0, 0.0);
    }

    update_e(&s);

    double scale = s.g->dt / (EPS0 * s.g->dx);

    for (size_t i = ibegin; i < iend; i++) {
        for (size_t j = jbegin; j < jend; j++) {
            double ex_aligned;
            double ex_cross;
            double ey_aligned;
            double ey_cross;

            seam_expected_curls(
                i, j,
                ax, ay, bx, by,
                &ex_aligned,
                &ex_cross,
                &ey_aligned,
                &ey_cross);

            double expected_ex = ex_old + scale * (
                                     mxx * ex_aligned + mxy * ex_cross);

            double expected_ey = ey_old + scale * (
                                     myy * ey_aligned + mxy * ey_cross);

            size_t idx = fdtd_index(s.g, i, j, 0);
            size_t before = failures;

            CLOSE(s.g->ex[idx], expected_ex, 1e-10, 1e-12);
            CLOSE(s.g->ey[idx], expected_ey, 1e-10, 1e-12);

            if (failures != before) {
                fprintf(
                    stderr,
                    "  seam context: axes=%u cpml=%d "
                    "mxy=%+.3f i=%zu j=%zu\n",
                    periodic_axes, with_cpml, mxy, i, j);
            }
        }
    }

cleanup:
    xt_destroy(&s);
}

static void test_periodic_tensor_seam_suite(void) {
    const double couplings[2] = {0.12, -0.12};

    for (size_t q = 0; q < 2; q++) {
        double mxy = couplings[q];

        /* Periodic x/z; nonperiodic y. */
        test_periodic_tensor_seams(5u, 0, mxy);
        test_periodic_tensor_seams(5u, 1, mxy);

        /* Periodic y/z; nonperiodic x. */
        test_periodic_tensor_seams(6u, 0, mxy);
        test_periodic_tensor_seams(6u, 1, mxy);

        /* Periodic x/y/z, including simultaneous seam crossings. */
        test_periodic_tensor_seams(7u, 0, mxy);
    }
}

static void test_uniform_tensor_evolution(
    unsigned mode_x,
    unsigned mode_y,
    double mxy)
{
    xt_sim s = {0};

    const size_t nx = 24;
    const size_t ny = 24;
    const size_t steps = 2048;

    const double mxx = 0.65;
    const double myy = 0.40;

    const double ax = 2.0 * XT_PI * (double)mode_x / (double)nx;
    const double ay = 2.0 * XT_PI * (double)mode_y / (double)ny;

    const double impedance = sqrt(MU0 / EPS0);

    /*
     * Nonzero initial electric and magnetic amplitudes exercise
     * the complete coupled recurrence from the first step.
     */
    double complex ex_ref = 0.13 + 0.07 * I;
    double complex ey_ref = -0.17 + 0.11 * I;
    double complex hz_ref = (0.9 - 0.2 * I) / impedance;

    double initial_size =
        cabs(ex_ref) + cabs(ey_ref) + impedance * cabs(hz_ref);

    double peak_size = initial_size;
    double max_normalized_error = 0.0;

    fprintf(
        stdout,
        "uniform tensor evolution: mode=(%u,%u) "
        "mxy=%+.3f steps=%zu\n",
        mode_x, mode_y, mxy, steps);
    fflush(stdout);

    REQUIRE(xt_create(
        &s, nx, ny, 1, 0.001, 7u, 0u) == 0);

    for (size_t i = 0; i < nx; i++) {
        for (size_t j = 0; j < ny; j++) {
            size_t idx = fdtd_index(s.g, i, j, 0);
            double complex phase =
                cexp(I * (ax * (double)i + ay * (double)j));

            s.g->ex[idx] = creal(ex_ref * phase);
            s.g->ey[idx] = creal(ey_ref * phase);
            s.g->hz[idx] = creal(hz_ref * phase);

            REQUIRE(materials_set_component_tensor(
                s.m, s.g, FIELD_EX, i, j, 0,
                mxx, mxy, myy) == 0);

            REQUIRE(materials_set_component_tensor(
                s.m, s.g, FIELD_EY, i, j, 0,
                myy, mxy, mxx) == 0);
        }
    }

    /*
     * Forward differences used by the magnetic update.
     */
    double complex forward_x = cexp(I * ax) - 1.0;
    double complex forward_y = cexp(I * ay) - 1.0;

    /*
     * Backward aligned differences used by the electric update.
     */
    double complex backward_x = 1.0 - cexp(-I * ax);
    double complex backward_y = 1.0 - cexp(-I * ay);

    /*
     * Ex cross curl: minus centered x difference, averaged in y.
     * Ey cross curl: plus centered y difference, averaged in x.
     */
    double complex cross_ex =
        -I * sin(ax) * 0.5 * (1.0 + cexp(-I * ay));

    double complex cross_ey =
         I * sin(ay) * 0.5 * (1.0 + cexp(-I * ax));

    double hscale = s.g->dt / (MU0 * s.g->dx);
    double escale = s.g->dt / (EPS0 * s.g->dx);

    for (size_t step = 1; step <= steps; step++) {
        /*
         * Independent modal recurrence, in the same leapfrog
         * update order as the grid: H first, then E.
         */
        hz_ref += hscale * (
            forward_y * ex_ref - forward_x * ey_ref);

        ex_ref += escale * (
            mxx * backward_y + mxy * cross_ex) * hz_ref;

        ey_ref += escale * (
            -myy * backward_x + mxy * cross_ey) * hz_ref;

        fdtd_update_h(s.g);
        update_e(&s);

        double modal_size =
            cabs(ex_ref) +
            cabs(ey_ref) +
            impedance * cabs(hz_ref);

        CHECK(isfinite(modal_size));

        if (!isfinite(modal_size)) {
            fprintf(
                stderr,
                "  nonfinite modal reference at step=%zu\n",
                step);
            goto cleanup;
        }

        if (modal_size > peak_size) {
            peak_size = modal_size;
        }

        /*
         * Generous finite-run growth guard, not an energy-
         * conservation assertion or a general stability theorem.
         */
        CHECK(modal_size <= 10.0 * initial_size);

        if (modal_size > 10.0 * initial_size) {
            fprintf(
                stderr,
                "  modal growth guard failed at step=%zu\n",
                step);
            goto cleanup;
        }

        if (step != 1 && step % 32 != 0 && step != steps) {
            continue;
        }

        size_t before = failures;

        for (size_t i = 0; i < nx; i++) {
            for (size_t j = 0; j < ny; j++) {
                size_t idx = fdtd_index(s.g, i, j, 0);
                double complex phase =
                    cexp(I * (ax * (double)i + ay * (double)j));

                double expected_ex = creal(ex_ref * phase);
                double expected_ey = creal(ey_ref * phase);
                double expected_hz =
                    impedance * creal(hz_ref * phase);

                double actual_hz = impedance * s.g->hz[idx];

                CLOSE(
                    s.g->ex[idx], expected_ex,
                    1e-8, 1e-10);

                CLOSE(
                    s.g->ey[idx], expected_ey,
                    1e-8, 1e-10);

                CLOSE(
                    actual_hz, expected_hz,
                    1e-8, 1e-10);

                /* The initially absent polarization stays absent. */
                CLOSE(s.g->ez[idx], 0.0, 0.0, 1e-12);
                CLOSE(impedance * s.g->hx[idx], 0.0, 0.0, 1e-12);
                CLOSE(impedance * s.g->hy[idx], 0.0, 0.0, 1e-12);

                double err =
                    fabs(s.g->ex[idx] - expected_ex);

                double ey_err =
                    fabs(s.g->ey[idx] - expected_ey);

                double hz_err =
                    fabs(actual_hz - expected_hz);

                if (ey_err > err) err = ey_err;
                if (hz_err > err) err = hz_err;

                err /= initial_size;

                if (err > max_normalized_error) {
                    max_normalized_error = err;
                }
            }
        }

        if (failures != before) {
            fprintf(
                stderr,
                "  evolution context: mode=(%u,%u) "
                "mxy=%+.3f step=%zu\n",
                mode_x, mode_y, mxy, step);
            goto cleanup;
        }
    }

    fprintf(
        stdout,
        "  modal peak/initial=%.6g "
        "max normalized field error=%.6g\n",
        peak_size / initial_size,
        max_normalized_error);

cleanup:
    xt_destroy(&s);
}

static void test_uniform_tensor_evolution_suite(void)
{
    /* Longer wavelength, both off-diagonal signs. */
    test_uniform_tensor_evolution(2u, 3u,  0.12);
    test_uniform_tensor_evolution(2u, 3u, -0.12);

    /* Shorter wavelength, away from the exact Nyquist null. */
    test_uniform_tensor_evolution(9u, 10u,  0.12);
    test_uniform_tensor_evolution(9u, 10u, -0.12);
}

static void variable_zero_fields(fdtd_grid_t *g)
{
    size_t bytes = grid_storage(g) * sizeof(double);

    memset(g->ex, 0, bytes);
    memset(g->ey, 0, bytes);
    memset(g->ez, 0, bytes);
    memset(g->hx, 0, bytes);
    memset(g->hy, 0, bytes);
    memset(g->hz, 0, bytes);
}

/*
 * pattern:
 *   0 = uniform control
 *   1 = x-directed two-region slab
 *   2 = periodically repeated tilted bands
 *   3 = cell-scale checkerboard
 *
 * Evaluate the region at each component's physical Yee location.
 */
static int variable_region(
    unsigned pattern,
    double x,
    double y,
    size_t n)
{
    switch (pattern) {
    case 0:
        return 0;

    case 1:
        return x >= 0.5 * (double)n;

    case 2:
        return fmod(x + 2.0 * y + 0.37, (double)n)
            >= 0.5 * (double)n;

    default:
        return ((size_t)floor(x) + (size_t)floor(y)) % 2u;
    }
}

static int variable_assign(
    xt_sim *s,
    unsigned pattern,
    double coupling_sign)
{
    for (size_t i = 0; i < s->g->nx; i++) {
        for (size_t j = 0; j < s->g->ny; j++) {
            for (int component = 0; component < 2; component++) {
                int is_ex = component == 0;

                double x = (double)i + (is_ex ? 0.5 : 0.0);
                double y = (double)j + (is_ex ? 0.0 : 0.5);

                int region = variable_region(
                    pattern, x, y, s->g->nx);

                double mxx = region ? 0.25 : 0.65;
                double myy = region ? 0.55 : 0.40;
                double mxy = coupling_sign *
                    (region ? -0.08 : 0.12);

                int rc;

                if (is_ex) {
                    rc = materials_set_component_tensor(
                        s->m, s->g, FIELD_EX, i, j, 0,
                        mxx, mxy, myy);
                } else {
                    rc = materials_set_component_tensor(
                        s->m, s->g, FIELD_EY, i, j, 0,
                        myy, mxy, mxx);
                }

                if (rc != 0) {
                    return -1;
                }
            }
        }
    }

    return 0;
}

/*
 * Probe the dimensionless Hz wave operator:
 *
 *     A = C Xi C^T
 *
 * For each Hz basis vector:
 *   E <- escale * Xi * curl(Hz)
 *   Hz <- Hz - hscale * curl(E)
 *
 * Therefore:
 *   A*Hz = -(Hz_after - Hz_before)/(hscale*escale).
 *
 * All operations on fields use the production update functions.
 */
static void test_variable_tensor_operator(
    unsigned pattern,
    double coupling_sign)
{
    xt_sim s = {0};
    double *a = NULL;

    const size_t n = 8;
    const size_t count = n * n;

    fprintf(
        stdout,
        "variable tensor operator: pattern=%u sign=%+.0f\n",
        pattern, coupling_sign);
    fflush(stdout);

    REQUIRE(xt_create(
        &s, n, n, 1, 0.001, 7u, 0u) == 0);

    REQUIRE(variable_assign(
        &s, pattern, coupling_sign) == 0);

    a = calloc(count * count, sizeof(*a));
    REQUIRE(a != NULL);

    double hscale = s.g->dt / (MU0 * s.g->dx);
    double escale = s.g->dt / (EPS0 * s.g->dx);
    double factor = hscale * escale;

    for (size_t column = 0; column < count; column++) {
        variable_zero_fields(s.g);

        size_t ci = column / n;
        size_t cj = column % n;

        s.g->hz[fdtd_index(s.g, ci, cj, 0)] = 1.0;

        update_e(&s);
        fdtd_update_h(s.g);

        for (size_t row = 0; row < count; row++) {
            size_t ri = row / n;
            size_t rj = row % n;

            double old_hz = row == column ? 1.0 : 0.0;
            double new_hz =
                s.g->hz[fdtd_index(s.g, ri, rj, 0)];

            a[row * count + column] =
                -(new_hz - old_hz) / factor;

            CHECK(isfinite(a[row * count + column]));
        }
    }

    double max_entry = 0.0;
    double max_asymmetry = 0.0;
    double max_row_sum = 0.0;

    for (size_t row = 0; row < count; row++) {
        double row_sum = 0.0;

        for (size_t column = 0; column < count; column++) {
            double value = a[row * count + column];

            double asymmetry = fabs(
                value - a[column * count + row]);

            if (fabs(value) > max_entry) {
                max_entry = fabs(value);
            }

            if (asymmetry > max_asymmetry) {
                max_asymmetry = asymmetry;
            }

            row_sum += value;
        }

        if (fabs(row_sum) > max_row_sum) {
            max_row_sum = fabs(row_sum);
        }
    }

    CHECK(max_entry > 0.0);
    CHECK(isfinite(max_asymmetry));
    CHECK(isfinite(max_row_sum));

    /*
     * A spatially constant Hz must have zero curl, even with
     * variable coefficients.
     */
    CHECK(max_row_sum <= 1e-11 * fmax(1.0, max_entry));

    /* Reciprocal assembly must be symmetric for every pattern. */
    CHECK(max_asymmetry <= 1e-11 * fmax(1.0, max_entry));

    fprintf(
        stdout,
        "  max|A|=%.9g relative asymmetry=%.9g "
        "constant-mode residual=%.9g\n",
        max_entry,
        max_asymmetry / fmax(max_entry, 1e-300),
        max_row_sum / fmax(max_entry, 1e-300));

cleanup:
    free(a);
    xt_destroy(&s);
}

/*
 * A finite-run stress test, not an energy-conservation test.
 * Monitor all six canonical field components.
 */
static void test_variable_tensor_evolution(
    unsigned pattern,
    double coupling_sign)
{
    xt_sim s = {0};

    const size_t n = 16;
    const size_t steps = 8192;
    const double impedance = sqrt(MU0 / EPS0);

    double initial_norm = 0.0;
    double peak_ratio = 1.0;

    fprintf(
        stdout,
        "variable tensor evolution: pattern=%u "
        "sign=%+.0f steps=%zu\n",
        pattern, coupling_sign, steps);
    fflush(stdout);

    REQUIRE(xt_create(
        &s, n, n, 1, 0.001, 7u, 0u) == 0);

    REQUIRE(variable_assign(
        &s, pattern, coupling_sign) == 0);

    /*
     * Deterministic, spatially varied initial data.
     * No sources are applied during evolution.
     */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            size_t idx = fdtd_index(s.g, i, j, 0);

            double value =
                sin(0.71 * (double)(i + 1) * (double)(j + 2)) +
                0.3 * cos(1.13 * (double)(2 * i + j + 1));

            s.g->hz[idx] = value / impedance;
            initial_norm += value * value;
        }
    }

    REQUIRE(isfinite(initial_norm) && initial_norm > 0.0);

    for (size_t step = 1; step <= steps; step++) {
        fdtd_update_h(s.g);
        update_e(&s);

        if (step != 1 && step % 32 != 0 && step != steps) {
            continue;
        }

        double norm = 0.0;
        double unwanted = 0.0;

        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                size_t idx = fdtd_index(s.g, i, j, 0);

                double ex = s.g->ex[idx];
                double ey = s.g->ey[idx];
                double ez = s.g->ez[idx];

                double hx = impedance * s.g->hx[idx];
                double hy = impedance * s.g->hy[idx];
                double hz = impedance * s.g->hz[idx];

                norm +=
                    ex * ex + ey * ey + ez * ez +
                    hx * hx + hy * hy + hz * hz;

                unwanted += ez * ez + hx * hx + hy * hy;
            }
        }

        double ratio = norm / initial_norm;

        CHECK(isfinite(ratio));

        if (!isfinite(ratio)) {
            fprintf(
                stderr,
                "  nonfinite fields at step=%zu\n", step);
            goto cleanup;
        }

        if (ratio > peak_ratio) {
            peak_ratio = ratio;
        }

        /*
         * Explicit regression guard:
         * 100x squared-field norm = 10x aggregate field amplitude.
         * This is not a derived physical energy bound.
         */
        CHECK(ratio <= 100.0);

        CHECK(unwanted <=
            1e-20 * fmax(initial_norm, norm));

        if (ratio > 100.0) {
            fprintf(
                stderr,
                "  growth guard failed: step=%zu "
                "norm/initial=%.9g\n",
                step, ratio);
            goto cleanup;
        }
    }

    fprintf(
        stdout,
        "  peak squared-field norm/initial=%.9g\n",
        peak_ratio);

cleanup:
    xt_destroy(&s);
}

static void test_variable_tensor_suite(void)
{
    /* Uniform off-diagonal control. */
    test_variable_tensor_operator(0u, 1.0);
    test_variable_tensor_evolution(0u, 1.0);

    /* Variable diagonal-only control. */
    test_variable_tensor_operator(1u, 0.0);
    test_variable_tensor_evolution(1u, 0.0);

    /* Slab, tilted bands, and checkerboard; both coupling signs. */
    for (unsigned pattern = 1u; pattern <= 3u; pattern++) {
        for (int sign = -1; sign <= 1; sign += 2) {
            test_variable_tensor_operator(pattern, (double)sign);
            test_variable_tensor_evolution(pattern, (double)sign);
        }
    }
}


static void test_isolated_tensor_edge(void)
{
    xt_sim s = {0};
    const size_t i = 10, j = 11;
    const double mxx = 0.65, mxy = 0.12, myy = 0.40;
    REQUIRE(create_fixture(&s, 0) == 0);
    for (size_t u = 0; u < s.g->nx; u++) {
        for (size_t v = 0; v < s.g->ny; v++) {
            s.g->hz[fdtd_index(s.g, u, v, 0)] =
                0.7 * (double)u - 0.2 * (double)v;
        }
    }
    REQUIRE(materials_set_component_tensor(
        s.m, s.g, FIELD_EX, i, j, 0, mxx, mxy, myy) == 0);
    update_e(&s);
    double scale = s.g->dt / (EPS0 * s.g->dx);
    CLOSE(s.g->ex[fdtd_index(s.g, i, j, 0)],
          scale * mxx * (-0.2), 1e-10, 1e-12);
    CLOSE(s.g->ey[fdtd_index(s.g, i, j, 0)],
          -scale * 0.7, 1e-10, 1e-12);
    fprintf(stdout, "isolated tensor edge: reciprocal neighbors absent\n");
cleanup:
    xt_destroy(&s);
}

int main(void) {
    test_api();

    test_analytic_update(0, 0);
    test_analytic_update(0, 1);
    test_analytic_update(1, 0);
    test_analytic_update(1, 1);

    test_diagonal_invariance(0);
    test_diagonal_invariance(1);

    test_periodic_tensor_seam_suite();
    test_uniform_tensor_evolution_suite();
    test_variable_tensor_suite();
    test_isolated_tensor_edge();

    printf(
        "materials_tensor_unit: %zu checks, %zu failures\n",
        checks, failures);

    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
