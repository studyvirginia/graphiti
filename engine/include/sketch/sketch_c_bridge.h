#ifndef SKETCH_C_BRIDGE_H
#define SKETCH_C_BRIDGE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int nodes;
    float shape_err;
    float r2;
    float dl_bits;
    const char* expr_x;
    const char* expr_y;
    int curve_count;
    const float* curve_x;
    const float* curve_y;
} SketchBridgeCandidate;

typedef void (*SketchBridgeCallback)(
    int gen,
    int elapsed_ms,
    const SketchBridgeCandidate* champion,
    const SketchBridgeCandidate* candidates,
    int candidate_count,
    void* user_data
);

typedef struct {
    int abi_version;
    const float* point_s2;
    const float* point_log_var;
    const float* point_nu;
    int point_count;
    float floor_clamp;
    float floor_stop;
    int floor_simplify_gens;
} SketchNoiseField;

#define SKETCH_NOISE_FIELD_ABI_VERSION 5

int sketch_solve_stroke_c(
    const float* pts_x,
    const float* pts_y,
    int n_pts,
    int pop,
    int gens,
    int max_length,
    SketchBridgeCallback callback,
    void* user_data
);

int sketch_solve_stroke_with_noise_c(
    const float* pts_x,
    const float* pts_y,
    int n_pts,
    int pop,
    int gens,
    int max_length,
    const SketchNoiseField* noise,
    SketchBridgeCallback callback,
    void* user_data
);

int sketch_gpu_pointseg_counters_c(long long* out4);
void sketch_gpu_pointseg_counters_reset_c(void);

void sketch_op_diag_print_c(void);

void sketch_nsafety_counters_print_c(void);

void sketch_op2_diag_print_c(void);

void sketch_e10_diag_print_c(void);

void sketch_diag_set_one_residual_c(int mode);
void sketch_diag_set_fp_seed_c(long long seed);

int sketch_diag_multitree_fingerprint_c(
    const char* id,
    const float* pts_x,
    const float* pts_y,
    int n_pts,
    int pop,
    int gens,
    int n_threads,
    int gpu_mode,
    char* out,
    int out_cap
);

void sketch_diag_set_cis_exact_c(int mode);
void sketch_diag_set_sf2_c(int mask);
int sketch_diag_multitree_rulers_c(
    const char* id,
    const float* pts_x,
    const float* pts_y,
    int n_pts,
    int pop,
    int gens,
    int n_threads,
    char* out,
    int out_cap
);

int sketch_sigma_v2_cells_c(
    const float* x,
    const float* y,
    int n,
    double* lnE_out,
    double* sd_out,
    int cap
);

#ifdef __cplusplus
}
#endif

#endif
