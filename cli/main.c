#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <complex.h>
#include <signal.h>
#include <sys/stat.h>
#include <errno.h>

#include "config.h"
#include "presets.h"
#include "settings.h"
#include "../engine/core/rf_signal.h"
#include "../engine/core/chain.h"

#if defined(_WIN32) && !defined(__CYGWIN__)
#include <direct.h>
#define MKDIR(path) _mkdir(path)
#else
#define MKDIR(path) mkdir(path, 0755)
#endif

#define MAX_STAGES 16
#define CONFIG_DIR "configs"
#define CONFIG_PATH "configs/config.json"
#define LOCK_PATH "configs/.rfsyn.lock"
#define STOP_PATH "configs/.rfsyn.stop"

static volatile sig_atomic_t g_stop_requested = 0;

static void handle_stop_signal(int sig) {
    (void)sig;
    g_stop_requested = 1;
}

static int ensure_directory(const char *path) {
    struct stat st;

    if (!path || !*path) return -1;

    if (MKDIR(path) == 0) return 0;

    if (errno != EEXIST || stat(path, &st) != 0) return -1;

#if defined(_WIN32) && !defined(__CYGWIN__)
    return (st.st_mode & _S_IFMT) == _S_IFDIR ? 0 : -1;
#else
    return S_ISDIR(st.st_mode) ? 0 : -1;
#endif
}

static signal_t *make_constant_signal(size_t n_samples) {
    signal_t *sig = signal_create(n_samples, 1e6, 915e6);
    if (!sig) return NULL;

    for (size_t i = 0; i < n_samples; i++) {
        sig->samples[i] = 1.0f + 0.0f * I;
    }

    return sig;
}

static int write_signal(const signal_t *sig, const char *path) {
    if (!sig || !path) return -1;
    if (sig->n_samples > 0 && !sig->samples) return -1;

    FILE *f = fopen(path, "wb");
    if (!f) return -1;

    size_t written = 0;
    if (sig->n_samples > 0) {
        written = fwrite(
            sig->samples,
            sizeof(*sig->samples),
            sig->n_samples,
            f
        );
    }

    int close_result = fclose(f);
    return (written == sig->n_samples && close_result == 0)
        ? 0 : -1;
}

static int run_is_active(void) {
    struct stat st;
    return stat(LOCK_PATH, &st) == 0;
}

static int create_lock(void) {
    FILE *f = fopen(LOCK_PATH, "wx");
    if (!f) return -1;

    if (fclose(f) != 0) {
        remove(LOCK_PATH);
        return -1;
    }

    return 0;
}

static void remove_lock(void) {
    remove(LOCK_PATH);
}

static int stop_file_present(void) {
    struct stat st;
    return stat(STOP_PATH, &st) == 0;
}

static void remove_stop_file(void) {
    remove(STOP_PATH);
}

static int cmd_init(void) {
    struct stat st;

    if (stat(CONFIG_PATH, &st) == 0) {
        fprintf(
            stderr,
            "%s already exists -- use "
            "`rfsyn config preset realistic` to reset values\n",
            CONFIG_PATH
        );
        return 1;
    }

    if (ensure_directory(CONFIG_DIR) != 0) {
        fprintf(stderr, "could not create %s/\n", CONFIG_DIR);
        return 1;
    }

    config_t *cfg = preset_build("realistic");
    if (!cfg) {
        fprintf(stderr, "failed to build default config\n");
        return 1;
    }

    const char *out_dir = config_get(cfg, "job", "output_dir");
    if (!out_dir) out_dir = PRESET_REALISTIC_JOB_OUTPUTDIR;

    if (ensure_directory(out_dir) != 0) {
        fprintf(stderr, "could not create output directory '%s'\n", out_dir);
        config_destroy(cfg);
        return 1;
    }

    if (config_save(cfg, CONFIG_PATH) != 0) {
        fprintf(stderr, "failed to write %s\n", CONFIG_PATH);
        config_destroy(cfg);
        return 1;
    }

    printf("created %s/\n", CONFIG_DIR);
    printf("wrote %s with realistic defaults:\n", CONFIG_PATH);
    config_print(cfg, stdout);
    printf("output directory: %s/\n", out_dir);
    printf("run `rfsyn start` to generate\n");

    config_destroy(cfg);
    return 0;
}

static int cmd_config_view(void) {
    config_t *cfg = config_load(CONFIG_PATH);
    if (!cfg) {
        fprintf(
            stderr,
            "could not load %s - check the file or run `rfsyn init`\n",
            CONFIG_PATH
        );
        return 1;
    }

    config_print(cfg, stdout);
    config_destroy(cfg);
    return 0;
}

static int cmd_config_preset(const char *name) {
    config_t *cfg = preset_build(name);
    if (!cfg) {
        fprintf(stderr, "unknown preset '%s' - available: realistic\n", name);
        return 1;
    }

    if (ensure_directory(CONFIG_DIR) != 0 ||
        config_save(cfg, CONFIG_PATH) != 0) {
        fprintf(stderr, "failed to write %s\n", CONFIG_PATH);
        config_destroy(cfg);
        return 1;
    }

    printf("applied preset '%s' to %s\n", name, CONFIG_PATH);
    config_destroy(cfg);
    return 0;
}

static int try_parse_bool(const char *s, int *out) {
    if (!strcmp(s, "true")) {
        *out = 1;
        return 1;
    }

    if (!strcmp(s, "false")) {
        *out = 0;
        return 1;
    }

    return 0;
}

static int cmd_config_set(const char *dotted_key, const char *value) {
    char section[128];
    char key[128];

    const char *dot = strchr(dotted_key, '.');
    if (!dot) {
        fprintf(stderr, "key must use section.key form\n");
        return 1;
    }

    size_t section_len = (size_t)(dot - dotted_key);
    size_t key_len = strlen(dot + 1);

    if (section_len == 0 || section_len >= sizeof(section) ||
        key_len == 0 || key_len >= sizeof(key)) {
        fprintf(stderr, "section/key is empty or too long\n");
        return 1;
    }

    memcpy(section, dotted_key, section_len);
    section[section_len] = '\0';
    memcpy(key, dot + 1, key_len + 1);

    config_t *cfg = config_load(CONFIG_PATH);
    if (!cfg) {
        fprintf(stderr, "could not load %s\n", CONFIG_PATH);
        return 1;
    }

    int bool_value;
    char *endptr = NULL;

    errno = 0;
    double double_value = strtod(value, &endptr);

    if (try_parse_bool(value, &bool_value)) {
        config_set_bool(cfg, section, key, bool_value);
    } else if (endptr != value && *endptr == '\0') {
        if (errno == ERANGE ||
            double_value != double_value ||
            double_value > 1.7976931348623157e308 ||
            double_value < -1.7976931348623157e308) {
            fprintf(stderr, "numeric value must be finite and in range\n");
            config_destroy(cfg);
            return 1;
        }

        config_set_double(cfg, section, key, double_value);
    } else {
        config_set_string(cfg, section, key, value);
    }

    int result = config_save(cfg, CONFIG_PATH);
    config_destroy(cfg);

    if (result != 0) {
        fprintf(stderr, "failed to write %s\n", CONFIG_PATH);
        return 1;
    }

    printf("set %s.%s = %s\n", section, key, value);
    return 0;
}

static int cmd_start(void) {
    config_t *cfg = config_load(CONFIG_PATH);
    if (!cfg) {
        fprintf(stderr, "could not load %s -- run `rfsyn init`\n", CONFIG_PATH);
        return 1;
    }

    long count = config_get_long(
        cfg, "job", "count", PRESET_REALISTIC_JOB_COUNT
    );

    long n_samples = config_get_long(
        cfg, "job", "n_samples", PRESET_REALISTIC_JOB_NSAMPLES
    );

    int fdtd_enabled = config_get_bool(
        cfg, "fdtd", "enabled", PRESET_FDTD_ENABLED
    );

    const char *out_dir = config_get(cfg, "job", "output_dir");
    if (!out_dir) out_dir = PRESET_REALISTIC_JOB_OUTPUTDIR;

    if (count <= 0) {
        fprintf(stderr, "job.count must be positive, got %ld\n", count);
        config_destroy(cfg);
        return 1;
    }

    if (!fdtd_enabled &&
        (n_samples <= 0 ||
         (uintmax_t)n_samples > SIZE_MAX / sizeof(float complex))) {
        fprintf(stderr, "job.n_samples must be positive and allocatable\n");
        config_destroy(cfg);
        return 1;
    }

    if (ensure_directory(CONFIG_DIR) != 0 ||
        ensure_directory(out_dir) != 0) {
        fprintf(stderr, "could not create config/output directories\n");
        config_destroy(cfg);
        return 1;
    }

    if (create_lock() != 0) {
        fprintf(
            stderr,
            "could not acquire %s - a job may already be running; "
            "if no job is running, check for a stale lock\n",
            LOCK_PATH
        );
        config_destroy(cfg);
        return 1;
    }

    g_stop_requested = 0;
    remove_stop_file();

    signal(SIGINT, handle_stop_signal);
    signal(SIGTERM, handle_stop_signal);

    printf(
        "generating %ld example(s) into '%s/' (%s)\n",
        count,
        out_dir,
        fdtd_enabled ? "FDTD simulation" : "constant-signal debug mode"
    );

    time_t last_report = time(NULL);
    long failures = 0;
    long attempted = 0;
    long written_ok = 0;

    for (long i = 0; i < count; i++) {
        if (g_stop_requested || stop_file_present()) {
            g_stop_requested = 1;
            break;
        }

        attempted++;

        transform_t stages[MAX_STAGES] = {0};
        size_t n_stages = 0;

        if (settings_build_chain(
                cfg, (uint64_t)i, stages, MAX_STAGES, &n_stages
            ) != 0) {
            fprintf(stderr, "could not build transform chain for example %ld\n", i);
            chain_free_stages(stages, n_stages);
            failures++;
            continue;
        }

        signal_t *sig = fdtd_enabled
            ? signal_create(0, 0.0, 0.0)
            : make_constant_signal((size_t)n_samples);

        if (!sig) {
            fprintf(stderr, "could not allocate signal for example %ld\n", i);
            chain_free_stages(stages, n_stages);
            failures++;
            continue;
        }

        if (chain_apply(stages, n_stages, sig) != 0) {
            fprintf(stderr, "chain failed on example %ld\n", i);
            failures++;
        } else {
            char path[1024];
            int path_len = snprintf(
                path, sizeof(path), "%s/example_%06ld.iq", out_dir, i
            );

            if (path_len < 0 || (size_t)path_len >= sizeof(path)) {
                fprintf(stderr, "output path too long for example %ld\n", i);
                failures++;
            } else if (write_signal(sig, path) != 0) {
                fprintf(stderr, "could not write '%s'\n", path);
                failures++;
            } else {
                written_ok++;
            }
        }

        signal_destroy(sig);
        chain_free_stages(stages, n_stages);

        time_t now = time(NULL);
        if (difftime(now, last_report) >= 1.0 || i == count - 1) {
            printf(
                "\r  %ld / %ld attempted; %ld written",
                attempted, count, written_ok
            );
            fflush(stdout);
            last_report = now;
        }
    }

    printf("\n");
    remove_lock();
    remove_stop_file();

    printf(
        "%s: %ld / %ld example(s) written to '%s/'\n",
        g_stop_requested ? "stopped early" : "done",
        written_ok,
        count,
        out_dir
    );

    if (failures > 0) {
        fprintf(stderr, "%ld example(s) failed\n", failures);
    }

    config_destroy(cfg);
    return failures > 0 ? 1 : 0;
}

static int cmd_end(void) {
    if (!run_is_active()) {
        fprintf(stderr, "no run in progress\n");
        return 1;
    }

    FILE *f = fopen(STOP_PATH, "w");
    if (!f) {
        fprintf(stderr, "could not write %s\n", STOP_PATH);
        return 1;
    }

    if (fclose(f) != 0) {
        fprintf(stderr, "could not finish writing %s\n", STOP_PATH);
        return 1;
    }

    printf("stop requested - the job will finish its current example and exit\n");
    return 0;
}

static void print_usage(const char *prog) {
    fprintf(
        stderr,
        "usage: %s <command> [args]\n\n"
        "commands:\n"
        "  init                         create default config\n"
        "  config set <key> <value>     set section.key\n"
        "  config preset <name>         apply a preset\n"
        "  config view                  show config\n"
        "  start                        generate signals\n"
        "  end                          request job stop\n"
        "  help                         show usage\n",
        prog
    );
}

int main(int argc, char **argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    const char *command = argv[1];

    if (!strcmp(command, "init")) return cmd_init();
    if (!strcmp(command, "start")) return cmd_start();
    if (!strcmp(command, "end")) return cmd_end();

    if (!strcmp(command, "help")) {
        print_usage(argv[0]);
        return 0;
    }

    if (!strcmp(command, "config")) {
        if (argc < 3) {
            print_usage(argv[0]);
            return 1;
        }

        const char *sub = argv[2];

        if (!strcmp(sub, "view") && argc == 3) {
            return cmd_config_view();
        }

        if (!strcmp(sub, "preset") && argc == 4) {
            return cmd_config_preset(argv[3]);
        }

        if (!strcmp(sub, "set") && argc == 5) {
            return cmd_config_set(argv[3], argv[4]);
        }

        fprintf(stderr, "invalid config command or argument count\n");
        print_usage(argv[0]);
        return 1;
    }

    fprintf(stderr, "unknown command '%s'\n", command);
    print_usage(argv[0]);
    return 1;
}