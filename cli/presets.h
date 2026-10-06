#ifndef ENGINE_CLI_PRESETS_H
#define ENGINE_CLI_PRESETS_H

#include "config.h"

#define PRESET_REALISTIC_JOB_COUNT      100
#define PRESET_REALISTIC_JOB_NSAMPLES   1024
#define PRESET_REALISTIC_JOB_OUTPUTDIR  "output"

#define PRESET_REALISTIC_AWGN_ENABLED   1
#define PRESET_REALISTIC_AWGN_SNRDB     15.0
#define PRESET_REALISTIC_AWGN_SEED      42

#define PRESET_FDTD_ENABLED     1
#define PRESET_FDTD_NX          50
#define PRESET_FDTD_NY          50
#define PRESET_FDTD_NZ          50
#define PRESET_FDTD_DX          1e-3
#define PRESET_FDTD_COURANT     0.99
#define PRESET_FDTD_NSTEPS      200
#define PRESET_FDTD_NPML        8
#define PRESET_FDTD_SRC_I       (-1)		// -1 = center
#define PRESET_FDTD_SRC_J       (-1)
#define PRESET_FDTD_SRC_K       (-1)
#define PRESET_FDTD_PROBE_I     (-1)		// -1 = center + (6,0,0
#define PRESET_FDTD_PROBE_J     (-1)
#define PRESET_FDTD_PROBE_K     (-1)
#define PRESET_FDTD_T0          40.0		// dt
#define PRESET_FDTD_TAU         8.0

config_t *preset_build(const char *name);

#endif