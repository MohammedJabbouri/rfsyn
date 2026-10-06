#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_WIN32) && !defined(__CYGWIN__)
#include <direct.h>
#define MKDIR(path) _mkdir(path)
#define CHDIR(path) _chdir(path)
#else
#include <sys/stat.h>
#include <unistd.h>
#define MKDIR(path) mkdir(path, 0755)
#define CHDIR(path) chdir(path)
#endif

static int run_cmd(const char *rfsyn, const char *args) {
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "\"%s\" %s", rfsyn, args);
    printf("+ %s\n", cmd);
    int rc = system(cmd);
    fflush(stdout);
    return rc;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <path-to-rfsyn-binary>\n", argv[0]);
        return 1;
    }
    const char *rfsyn = argv[1];

    MKDIR("cli_smoke_work");
    if (CHDIR("cli_smoke_work") != 0) {
        perror("chdir");
        return 1;
    }

    remove("configs/config.json");
    remove("configs/.rfsyn.lock");
    remove("configs/.rfsyn.stop");

    if (run_cmd(rfsyn, "init") != 0) { fprintf(stderr, "FAIL: init\n"); return 1; }
    if (run_cmd(rfsyn, "config set job.count 5") != 0) { fprintf(stderr, "FAIL: config set job.count\n"); return 1; }
    if (run_cmd(rfsyn, "config set job.n_samples 256") != 0) { fprintf(stderr, "FAIL: config set job.n_samples\n"); return 1; }
    if (run_cmd(rfsyn, "config view") != 0) { fprintf(stderr, "FAIL: config view\n"); return 1; }
    if (run_cmd(rfsyn, "start") != 0) { fprintf(stderr, "FAIL: start\n"); return 1; }

    for (int i = 0; i < 5; i++) {
        char path[128];
        snprintf(path, sizeof(path), "output/example_%06d.iq", i);
        FILE *f = fopen(path, "rb");
        if (!f) {
            fprintf(stderr, "FAIL: %s missing\n", path);
            return 1;
        }
        if (i == 0) {
            float re = 0.0f, im = 0.0f;
            size_t got_re = fread(&re, sizeof(float), 1, f);
            size_t got_im = fread(&im, sizeof(float), 1, f);
            if (got_re != 1 || got_im != 1) {
                fprintf(stderr, "FAIL: could not read first sample of %s\n", path);
                fclose(f);
                return 1;
            }
            if (re == 1.0f && im == 0.0f) {
                fprintf(stderr, "FAIL: first sample of %s is exactly 1+0j -- AWGN did not run\n", path);
                fclose(f);
                return 1;
            }
            printf("OK: %s first sample = %.6f%+.6fj (AWGN applied)\n", path, re, im);
        }
        fclose(f);
    }

    printf("PASS\n");
    return 0;
}