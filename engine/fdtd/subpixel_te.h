#ifndef RFSYN_SUBPIXEL_TE_H
#define RFSYN_SUBPIXEL_TE_H

typedef enum {
    SUBPIXEL_TE_OK = 0,
    SUBPIXEL_TE_INVALID_ARGUMENT = -1,
    SUBPIXEL_TE_NUMERIC_ERROR = -2
} subpixel_te_status_t;

typedef struct {
    double center_x;
    double center_y;
    double size_u;
    double size_v;
    double angle_rad;
} subpixel_te_rectangle_t;

int subpixel_te_rectangle_fraction(
    const subpixel_te_rectangle_t *rectangle,
    double sample_x,
    double sample_y,
    double dx,
    double dy,
    double *fraction_out
);

int subpixel_te_mix_epsilon(
    double fraction,
    double eps_background,
    double eps_inside,
    double *epsilon_r_out
);

int subpixel_te_rectangle_epsilon(
    const subpixel_te_rectangle_t *rectangle,
    double sample_x,
    double sample_y,
    double dx,
    double dy,
    double eps_background,
    double eps_inside,
    double *fraction_out,
    double *epsilon_r_out
);

#endif