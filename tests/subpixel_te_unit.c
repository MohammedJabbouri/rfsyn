#include "../engine/fdtd/subpixel_te.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stddef.h>

#define TEST_PI 3.1415926535897932384626433832795

static size_t checks = 0;
static size_t failures = 0;

#define CHECK(expression) do { \
    checks++; \
    if (!(expression)) { \
        fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #expression); \
        failures++; \
    } \
} while (0)

static void close_to(
    double actual,
    double expected,
    double tolerance,
    int line)
{
    checks++;

    if (!isfinite(actual) ||
        fabs(actual - expected) > tolerance) {
        fprintf(
            stderr,
            "FAIL line %d: actual %.17g, expected %.17g, tolerance %.3g\n",
            line,
            actual,
            expected,
            tolerance
        );

        failures++;
    }
}

#define CLOSE(actual, expected, tolerance) \
    close_to((actual), (expected), (tolerance), __LINE__)

static double fraction(
    const subpixel_te_rectangle_t *rectangle,
    double x,
    double y,
    double dx,
    double dy)
{
    double value = NAN;

    CHECK(
        subpixel_te_rectangle_fraction(
            rectangle,
            x,
            y,
            dx,
            dy,
            &value
        ) == SUBPIXEL_TE_OK
    );

    return value;
}

static void test_basic_geometry(void)
{
    subpixel_te_rectangle_t rectangle = {
        0.0, 0.0, 4.0, 4.0, 0.0
    };

    CLOSE(fraction(&rectangle, 0, 0, 1, 1), 1, 1e-14);
    CLOSE(fraction(&rectangle, 10, 0, 1, 1), 0, 1e-14);

    rectangle = (subpixel_te_rectangle_t){
        1, 0, 2, 4, 0
    };

    CLOSE(fraction(&rectangle, 0, 0, 1, 1), 0.5, 1e-14);

    rectangle.center_y = 1;
    rectangle.size_v = 2;

    CLOSE(fraction(&rectangle, 0, 0, 1, 1), 0.25, 1e-14);

    rectangle = (subpixel_te_rectangle_t){
        0, 0, 0.5, 0.25, 0
    };

    CLOSE(fraction(&rectangle, 0, 0, 1, 1), 0.125, 1e-14);

    rectangle.angle_rad = TEST_PI / 6;

    CLOSE(fraction(&rectangle, 0, 0, 1, 1), 0.125, 1e-14);

    rectangle = (subpixel_te_rectangle_t){
        0, 0, 1, 1, 0
    };

    CLOSE(fraction(&rectangle, 0, 0, 2, 1), 0.5, 1e-14);

    rectangle = (subpixel_te_rectangle_t){
        1, 0, 1, 2, 0
    };

    CLOSE(fraction(&rectangle, 0, 0, 1, 1), 0, 1e-14);

    rectangle = (subpixel_te_rectangle_t){
        0,
        0,
        sqrt(0.5),
        sqrt(0.5),
        TEST_PI / 4
    };

    CLOSE(fraction(&rectangle, 0, 0, 1, 1), 0.5, 2e-14);

    rectangle = (subpixel_te_rectangle_t){
        1000001,
        -1000000,
        2,
        4,
        0
    };

    CLOSE(
        fraction(
            &rectangle,
            1000000,
            -1000000,
            1,
            1
        ),
        0.5,
        1e-14
    );
}

static void test_rotated_cuts(void)
{
    const double angles[] = {
        0, 15, 30, 45, 60, 90, 120
    };

    size_t count =
        sizeof(angles) / sizeof(angles[0]);

    for (size_t i = 0; i < count; i++) {
        double angle =
            angles[i] * TEST_PI / 180.0;

        subpixel_te_rectangle_t rectangle = {
            cos(angle),
            sin(angle),
            2,
            10,
            angle
        };

        CLOSE(
            fraction(&rectangle, 0, 0, 1, 1),
            0.5,
            1e-13
        );
    }

    const double cuts[] = {
        -1, -0.5, 0, 0.5, 1
    };

    const double expected[] = {
        0, 0.125, 0.5, 0.875, 1
    };

    for (size_t i = 0; i < 5; i++) {
        double center_u =
            cuts[i] / sqrt(2.0) - 2;

        subpixel_te_rectangle_t rectangle = {
            center_u / sqrt(2.0),
            center_u / sqrt(2.0),
            4,
            10,
            TEST_PI / 4
        };

        CLOSE(
            fraction(&rectangle, 0, 0, 1, 1),
            expected[i],
            1e-13
        );
    }
}

static void test_translation_and_units(void)
{
    double previous = 1.0;

    for (int q = -4; q <= 4; q++) {
        double shift = q * 0.1;

        subpixel_te_rectangle_t rectangle = {
            1 + shift,
            0,
            2,
            4,
            0
        };

        double value =
            fraction(&rectangle, 0, 0, 1, 1);

        CLOSE(value, 0.5 - shift, 1e-14);
        CHECK(value <= previous + 1e-14);

        previous = value;
    }

    subpixel_te_rectangle_t rectangle = {
        0.2,
        -0.1,
        0.7,
        0.6,
        0.43
    };

    double original =
        fraction(&rectangle, 0, 0, 1, 1);

    const double scales[] = {
        1e-3,
        1e3
    };

    for (size_t q = 0; q < 2; q++) {
        double scale = scales[q];

        subpixel_te_rectangle_t scaled = {
            rectangle.center_x * scale,
            rectangle.center_y * scale,
            rectangle.size_u * scale,
            rectangle.size_v * scale,
            rectangle.angle_rad
        };

        CLOSE(
            fraction(
                &scaled,
                0,
                0,
                scale,
                scale
            ),
            original,
            1e-13
        );
    }

    rectangle.angle_rad += 2 * TEST_PI;

    CLOSE(
        fraction(&rectangle, 0, 0, 1, 1),
        original,
        1e-13
    );
}

static void test_area_partition(void)
{
    for (int q = 0; q < 4; q++) {
        subpixel_te_rectangle_t rectangle = {
            0.13,
            -0.17,
            0.8,
            1.2,
            q * TEST_PI / 12
        };

        double area = 0.0;

        for (int i = 0; i < 16; i++) {
            for (int j = 0; j < 16; j++) {
                double x =
                    -2 + (i + 0.5) * 0.25;

                double y =
                    -2 + (j + 0.5) * 0.25;

                area +=
                    fraction(
                        &rectangle,
                        x,
                        y,
                        0.25,
                        0.25
                    ) * 0.25 * 0.25;
            }
        }

        CLOSE(area, 0.96, 2e-12);
    }
}

static uint32_t random_state =
    UINT32_C(123456789);

static double random_unit(void)
{
    random_state =
        UINT32_C(1664525) * random_state
        + UINT32_C(1013904223);

    return
        (double)(random_state >> 8)
        / 16777216.0;
}

static void test_subdivision(void)
{
    for (int q = 0; q < 500; q++) {
        subpixel_te_rectangle_t rectangle = {
            random_unit() - 0.5,
            random_unit() - 0.5,
            0.01 + 1.5 * random_unit(),
            0.01 + 1.5 * random_unit(),
            2 * TEST_PI * random_unit()
        };

        double coarse =
            fraction(&rectangle, 0, 0, 1, 1);

        double fine = 0.0;

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                fine += 0.25 * fraction(
                    &rectangle,
                    i ? 0.25 : -0.25,
                    j ? 0.25 : -0.25,
                    0.5,
                    0.5
                );
            }
        }

        CHECK(coarse >= 0.0 && coarse <= 1.0);
        CLOSE(fine, coarse, 1e-12);
    }
}

static void test_epsilon(void)
{
    double epsilon = -1.0;

    CHECK(
        subpixel_te_mix_epsilon(
            0, 1, 4, &epsilon
        ) == SUBPIXEL_TE_OK
    );

    CLOSE(epsilon, 1, 0);

    CHECK(
        subpixel_te_mix_epsilon(
            1, 1, 4, &epsilon
        ) == SUBPIXEL_TE_OK
    );

    CLOSE(epsilon, 4, 0);

    CHECK(
        subpixel_te_mix_epsilon(
            0.5, 1, 4, &epsilon
        ) == SUBPIXEL_TE_OK
    );

    CLOSE(epsilon, 2.5, 1e-14);

    CHECK(
        subpixel_te_mix_epsilon(
            0.25, 4, 1, &epsilon
        ) == SUBPIXEL_TE_OK
    );

    CLOSE(epsilon, 3.25, 1e-14);

    CHECK(
        subpixel_te_mix_epsilon(
            0.37, 7, 7, &epsilon
        ) == SUBPIXEL_TE_OK
    );

    CLOSE(epsilon, 7, 1e-14);

    CHECK(
        subpixel_te_mix_epsilon(
            0.5, 1, DBL_MAX, &epsilon
        ) == SUBPIXEL_TE_OK
    );

    CHECK(isfinite(epsilon));
    CLOSE(epsilon / DBL_MAX, 0.5, 1e-14);

    subpixel_te_rectangle_t rectangle = {
        1, 0, 2, 4, 0
    };

    double value = -1.0;

    CHECK(
        subpixel_te_rectangle_epsilon(
            &rectangle,
            0,
            0,
            1,
            1,
            1,
            4,
            &value,
            &epsilon
        ) == SUBPIXEL_TE_OK
    );

    CLOSE(value, 0.5, 1e-14);
    CLOSE(epsilon, 2.5, 1e-14);
}

static void test_errors(void)
{
    subpixel_te_rectangle_t rectangle = {
        0, 0, 1, 1, 0
    };

    double value = 0.123;
    double epsilon = 9.0;

    CHECK(
        subpixel_te_rectangle_fraction(
            NULL, 0, 0, 1, 1, &value
        ) == SUBPIXEL_TE_INVALID_ARGUMENT
    );

    CHECK(
        subpixel_te_rectangle_fraction(
            &rectangle, 0, 0, 1, 1, NULL
        ) == SUBPIXEL_TE_INVALID_ARGUMENT
    );

    CHECK(
        subpixel_te_rectangle_fraction(
            &rectangle, NAN, 0, 1, 1, &value
        ) == SUBPIXEL_TE_INVALID_ARGUMENT
    );

    CHECK(
        subpixel_te_rectangle_fraction(
            &rectangle, 0, 0, 0, 1, &value
        ) == SUBPIXEL_TE_INVALID_ARGUMENT
    );

    CHECK(
        subpixel_te_rectangle_fraction(
            &rectangle, 0, 0, 1, INFINITY, &value
        ) == SUBPIXEL_TE_INVALID_ARGUMENT
    );

    rectangle.size_u = -1;

    CHECK(
        subpixel_te_rectangle_fraction(
            &rectangle, 0, 0, 1, 1, &value
        ) == SUBPIXEL_TE_INVALID_ARGUMENT
    );

    rectangle.size_u = 1;
    rectangle.angle_rad = NAN;

    CHECK(
        subpixel_te_rectangle_fraction(
            &rectangle, 0, 0, 1, 1, &value
        ) == SUBPIXEL_TE_INVALID_ARGUMENT
    );

    CLOSE(value, 0.123, 0);

    rectangle.angle_rad = 0;

    CHECK(
        subpixel_te_mix_epsilon(
            -0.01, 1, 4, &epsilon
        ) == SUBPIXEL_TE_INVALID_ARGUMENT
    );

    CHECK(
        subpixel_te_mix_epsilon(
            1.01, 1, 4, &epsilon
        ) == SUBPIXEL_TE_INVALID_ARGUMENT
    );

    CHECK(
        subpixel_te_mix_epsilon(
            NAN, 1, 4, &epsilon
        ) == SUBPIXEL_TE_INVALID_ARGUMENT
    );

    CHECK(
        subpixel_te_mix_epsilon(
            0.5, 0.5, 4, &epsilon
        ) == SUBPIXEL_TE_INVALID_ARGUMENT
    );

    CHECK(
        subpixel_te_mix_epsilon(
            0.5, 1, INFINITY, &epsilon
        ) == SUBPIXEL_TE_INVALID_ARGUMENT
    );

    CLOSE(epsilon, 9, 0);

    CHECK(
        subpixel_te_rectangle_epsilon(
            &rectangle,
            0,
            0,
            1,
            1,
            1,
            4,
            &value,
            &value
        ) == SUBPIXEL_TE_INVALID_ARGUMENT
    );

    CHECK(
        subpixel_te_rectangle_epsilon(
            &rectangle,
            0,
            0,
            1,
            1,
            1,
            0,
            &value,
            &epsilon
        ) == SUBPIXEL_TE_INVALID_ARGUMENT
    );

    CLOSE(value, 0.123, 0);
    CLOSE(epsilon, 9, 0);

    rectangle.center_x = -DBL_MAX;

    CHECK(
        subpixel_te_rectangle_fraction(
            &rectangle,
            DBL_MAX,
            0,
            1,
            1,
            &value
        ) == SUBPIXEL_TE_NUMERIC_ERROR
    );

    CLOSE(value, 0.123, 0);
}

int main(void)
{
    test_basic_geometry();
    test_rotated_cuts();
    test_translation_and_units();
    test_area_partition();
    test_subdivision();
    test_epsilon();
    test_errors();

    printf(
        "subpixel_te_unit: %zu checks, %zu failures\n",
        checks,
        failures
    );

    return failures ? 1 : 0;
}