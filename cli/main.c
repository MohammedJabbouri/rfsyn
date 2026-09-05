#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <complex.h>
#include <signal.h>
#include <sys/stat.h>

#include "config.h"
#include "presets.h"
#include "settings.h"
#include "../engine/core/signal.h"
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

static signal_t *make_placeholder_signal(size_t n_samples) {
    signal_t *sig = signal_create(n_samples, 1e6, 915e6);
    if (!sig) return NULL;
    for (size_t i = 0; i < n_samples; i++) {
        sig->samples[i] = 1.0f + 0.0f * I;
    }
    return sig;
}

static int write_signal(const signal_t *sig, const char *path) {
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    size_t written = fwrite(sig->samples, sizeof(float complex), sig->n_samples, f);
    fclose(f);
    return (written == sig->n_samples) ? 0 : -1;
}

static int run_is_active(void) {
    FILE *f = fopen(LOCK_PATH, "r");
    if (f) { fclose(f); return 1; }
    return 0;
}

static void create_lock(void) {
    FILE *f = fopen(LOCK_PATH, "w");
    if (f) fclose(f);
}

static void remove_lock(void) {
    remove(LOCK_PATH);
}

static int stop_file_present(void) {
    FILE *f = fopen(STOP_PATH, "r");
    if (f) { fclose(f); return 1; }
    return 0;
}

static void remove_stop_file(void) {
    remove(STOP_PATH);
}

static int cmd_init(void) {
    struct stat st;
    if (stat(CONFIG_PATH, &st) == 0) {
        fprintf(stderr, "%s already exists -- delete it first, or use `rfsyn config preset realistic` to reset values\n", CONFIG_PATH);
        return 1;
    }

    MKDIR(CONFIG_DIR);

    config_t *cfg = preset_build("realistic");
    if (!cfg) {
        fprintf(stderr, "failed to build default config\n");
        return 1;
    }

    if (config_save(cfg, CONFIG_PATH) != 0) {
        fprintf(stderr, "failed to write %s\n", CONFIG_PATH);
        config_destroy(cfg);
        return 1;
    }

    const char *out_dir = config_get(cfg, "job", "output_dir");
    if (!out_dir) out_dir = PRESET_REALISTIC_JOB_OUTPUTDIR;
    MKDIR(out_dir);

    printf("created %s/\n", CONFIG_DIR);
    printf("wrote %s with realistic defaults:\n", CONFIG_PATH);
    printf("    job.count = %d\n", PRESET_REALISTIC_JOB_COUNT);
    printf("    job.n_samples = %d\n", PRESET_REALISTIC_JOB_NSAMPLES);
    printf("    job.output_dir = %s\n", PRESET_REALISTIC_JOB_OUTPUTDIR);
    printf("    awgn.enabled = true\n");
    printf("    awgn.snr_db = %.1f\n", PRESET_REALISTIC_AWGN_SNRDB);
    printf("    awgn.seed = %d\n", PRESET_REALISTIC_AWGN_SEED);
    printf("created %s/\n", out_dir);
    printf("run `rfsyn start` to generate\n");

    config_destroy(cfg);
    return 0;
}

static int cmd_config_view(void) {
    config_t *cfg = config_load(CONFIG_PATH);
    if (!cfg) {
        fprintf(stderr, "no config found at %s -- run `rfsyn init` first\n", CONFIG_PATH);
        return 1;
    }
    config_print(cfg, stdout);
    config_destroy(cfg);
    return 0;
}

static int cmd_config_preset(const char *name) {
    config_t *cfg = preset_build(name);
    if (!cfg) {
        fprintf(stderr, "unknown preset '%s' -- available: realistic\n", name);
        return 1;
    }

    MKDIR(CONFIG_DIR);
    if (config_save(cfg, CONFIG_PATH) != 0) {
        fprintf(stderr, "failed to write %s\n", CONFIG_PATH);
        config_destroy(cfg);
        return 1;
    }

    printf("applied preset '%s' to %s\n", name, CONFIG_PATH);
    config_destroy(cfg);
    return 0;
}

static int try_parse_bool(const char *s, int *out) {
    if (!strcmp(s, "true")) { *out = 1; return 1; }
    if (!strcmp(s, "false")) { *out = 0; return 1; }
    return 0;
}

static int cmd_config_set(const char *dotted_key, const char *value) {
    char section[128];
    char key[128];
    const char *dot = strchr(dotted_key, '.');
    if (!dot) {
        fprintf(stderr, "key must be in 'section.key' form, e.g. awgn.snr_db\n");
        return 1;
    }
    size_t section_len = (size_t)(dot - dotted_key);
    if (section_len >= sizeof(section)) {
        fprintf(stderr, "section name too long\n");
        return 1;
    }
    memcpy(section, dotted_key, section_len);
    section[section_len] = '\0';
    strncpy(key, dot + 1, sizeof(key) - 1);
    key[sizeof(key) - 1] = '\0';

    config_t *cfg = config_load(CONFIG_PATH);
    if (!cfg) {
        fprintf(stderr, "no config found at %s -- run `rfsyn init` first\n", CONFIG_PATH);
        return 1;
    }

    int bool_value;
    char *endptr;
    double double_value = strtod(value, &endptr);

    if (try_parse_bool(value, &bool_value)) {
        config_set_bool(cfg, section, key, bool_value);
    } else if (endptr != value && *endptr == '\0') {
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
    if (run_is_active()) {
        fprintf(stderr, "a run already appears to be in progress - use `rfsyn end`, or delete %s if you're sure nothing is running\n", LOCK_PATH);
        return 1;
    }

    config_t *cfg = config_load(CONFIG_PATH);
    if (!cfg) {
        fprintf(stderr, "no config found at %s - run `rfsyn init` first\n", CONFIG_PATH);
        return 1;
    }

    long count = config_get_long(cfg, "job", "count", PRESET_REALISTIC_JOB_COUNT);
    long n_samples = config_get_long(cfg, "job", "n_samples", PRESET_REALISTIC_JOB_NSAMPLES);
    const char *out_dir = config_get(cfg, "job", "output_dir");
    if (!out_dir) out_dir = PRESET_REALISTIC_JOB_OUTPUTDIR;

    if (count <= 0) {
        fprintf(stderr, "job.count must be positive, got %ld\n", count);
        config_destroy(cfg);
        return 1;
    }

    MKDIR(out_dir);
    signal(SIGINT, handle_stop_signal);
    signal(SIGTERM, handle_stop_signal);
    create_lock();
    remove_stop_file();

    printf("generating %ld example(s) into '%s/'\n", count, out_dir);

    time_t last_report = time(NULL);
    int failures = 0;
    long completed = 0;

    for (long i = 0; i < count; i++) {
        if (g_stop_requested || stop_file_present()) {
            g_stop_requested = 1;
            break;
        }

        transform_t stages[MAX_STAGES];
        size_t n_stages;

        if (settings_build_chain(cfg, (uint64_t)i, stages, MAX_STAGES, &n_stages) != 0) {
            fprintf(stderr, "too many stages enabled for MAX_STAGES=%d\n", MAX_STAGES);
            chain_free_stages(stages, n_stages);
            failures++;
            continue;
        }

        signal_t *sig = make_placeholder_signal((size_t)n_samples);
        if (!sig) {
            fprintf(stderr, "out of memory building example %ld\n", i);
            chain_free_stages(stages, n_stages);
            failures++;
            continue;
        }

        if (chain_apply(stages, n_stages, sig) != 0) {
            fprintf(stderr, "chain failed on example %ld\n", i);
            failures++;
        } else {
            char path[512];
            snprintf(path, sizeof(path), "%s/example_%06ld.iq", out_dir, i);
            if (write_signal(sig, path) != 0) {
                fprintf(stderr, "could not write '%s'\n", path);
                failures++;
            }
        }

        signal_destroy(sig);
        chain_free_stages(stages, n_stages);
        completed++;

        time_t now = time(NULL);
        if (difftime(now, last_report) >= 1.0 || i == count - 1) {
            printf("\r  %ld / %ld", i + 1, count);
            fflush(stdout);
            last_report = now;
        }
    }
    printf("\n");

    remove_lock();
    remove_stop_file();
    config_destroy(cfg);

    if (g_stop_requested) {
        printf("stopped early: %ld / %ld example(s) written to '%s/'\n", completed, count, out_dir);
        return 0;
    }
    if (failures > 0) {
        fprintf(stderr, "done with %d failure(s)\n", failures);
        return 1;
    }
    printf("done: %ld example(s) written to '%s/'\n", count, out_dir);
    return 0;
}

static int cmd_end(void) {
    if (!run_is_active()) {
        fprintf(stderr, "no run in progress\n");
        return 1;
    }

    FILE *f = fopen(STOP_PATH, "w");
    if (!f) {
        fprintf(stderr, "couldn't write %s\n", STOP_PATH);
        return 1;
    }
    fclose(f);

    printf("stop requested - the running job will finish its current example and exit\n");
    return 0;
}

static void print_usage(const char *prog) {
    fprintf(stderr,
        "usage: %s <command> [args]\n\n"
        "commands:\n"
        "  init                         create configs/config.json with realistic defaults\n"
        "  config set <key> <value>     set a config value, e.g. awgn.snr_db 12\n"
        "  config preset <name>         apply a named preset (currently: realistic)\n"
        "  config view                  print the current config\n"
        "  start                        run a generation job from configs/config.json\n"
        "  end                          stop a run started with `start`\n"
        "  help                         show this message\n",
        prog);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    const char *command = argv[1];

    if (!strcmp(command, "init")) return cmd_init();

    if (!strcmp(command, "config")) {
        if (argc < 3) {
            fprintf(stderr, "usage: %s config <set|preset|view> ...\n", argv[0]);
            return 1;
        }
        const char *sub = argv[2];
        if (!strcmp(sub, "set")) {
            if (argc != 5) {
                fprintf(stderr, "usage: %s config set <key> <value>\n", argv[0]);
                return 1;
            }
            return cmd_config_set(argv[3], argv[4]);
        }
        if (!strcmp(sub, "preset")) {
            if (argc != 4) {
                fprintf(stderr, "usage: %s config preset <name>\n", argv[0]);
                return 1;
            }
            return cmd_config_preset(argv[3]);
        }
        if (!strcmp(sub, "view")) return cmd_config_view();

        fprintf(stderr, "unknown config subcommand '%s'\n", sub);
        return 1;
    }

    if (!strcmp(command, "start")) return cmd_start();
    if (!strcmp(command, "end")) return cmd_end();

    if (!strcmp(command, "help") || !strcmp(command, "--help") || !strcmp(command, "-h")) {
        print_usage(argv[0]);
        return 0;
    }

    fprintf(stderr, "unknown command '%s'\n", command);
    print_usage(argv[0]);
    return 1;
}