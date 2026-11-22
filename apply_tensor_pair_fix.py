#!/usr/bin/env python3
"""Apply the reciprocal tensor-coupling candidate locally. No Git operations."""
from pathlib import Path
import re
import shutil
import sys
import time

root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
mp = root / "engine/fdtd/materials.c"
tp = root / "tests/materials_tensor_unit.c"
m = mp.read_text()
t = tp.read_text()
if "static double tensor_pair_weight(" in m:
    raise SystemExit("Pair patch already present; no files changed.")

def once(text, pattern, replacement, label):
    result, count = re.subn(pattern, lambda match: replacement, text, flags=re.S)
    if count != 1:
        raise SystemExit(f"Expected exactly one {label}; found {count}. No files changed.")
    return result

CBLOCK = '/*\n * Reciprocal tensor coupling.\n *\n * Cross-couplings exist only between two tensor-enabled edges.\n * Each Ex/Ey pair uses the same coefficient in both directions.\n * The diagonal entries remain component-local m11 values.\n */\nstatic double tensor_pair_weight(\n    const material_tensor_t *ex,\n    const material_tensor_t *ey)\n{\n    double rx = (ex->m12 / sqrt(ex->m11)) / sqrt(ex->m22);\n    double ry = (ey->m12 / sqrt(ey->m11)) / sqrt(ey->m22);\n    double rho = 0.5 * rx + 0.5 * ry;\n\n    return (0.25 * rho * sqrt(ex->m11)) * sqrt(ey->m11);\n}\n\nstatic inline void update_tensor_edge(\n    fdtd_grid_t *g,\n    const materials_t *mat,\n    const material_tensor_t *t,\n    int component,\n    size_t i,\n    size_t j,\n    size_t k)\n{\n    double value;\n\n    if (component == 0) {\n        size_t ip = material_next(i, g->nx);\n        size_t jm = material_previous(j, g->ny);\n        size_t ni[4] = { i, ip, i, ip };\n        size_t nj[4] = { j, j, jm, jm };\n\n        value = t->m11 * (HZ(i, j, k) - HZ(i, jm, k));\n\n        for (size_t q = 0; q < 4; q++) {\n            size_t u = ni[q];\n            size_t v = nj[q];\n            const material_tensor_t *other = tensor_lookup(\n                mat, fdtd_index(g, u, v, k), 1);\n\n            if (!other) continue;\n\n            size_t um = material_previous(u, g->nx);\n            double curl_y = -(HZ(u, v, k) - HZ(um, v, k));\n            value += tensor_pair_weight(t, other) * curl_y;\n        }\n\n        EX(i, j, k) += t->scale * value;\n    } else {\n        size_t im = material_previous(i, g->nx);\n        size_t jp = material_next(j, g->ny);\n        size_t ni[4] = { i, im, i, im };\n        size_t nj[4] = { j, j, jp, jp };\n\n        value = -t->m11 * (HZ(i, j, k) - HZ(im, j, k));\n\n        for (size_t q = 0; q < 4; q++) {\n            size_t u = ni[q];\n            size_t v = nj[q];\n            const material_tensor_t *other = tensor_lookup(\n                mat, fdtd_index(g, u, v, k), 0);\n\n            if (!other) continue;\n\n            size_t vm = material_previous(v, g->ny);\n            double curl_x = HZ(u, v, k) - HZ(u, vm, k);\n            value += tensor_pair_weight(other, t) * curl_x;\n        }\n\n        EY(i, j, k) += t->scale * value;\n    }\n}\n\n'

m = once(m,
    r"static inline void update_tensor_edge\(.*?(?=/\*\s*\* A NULL pml selects ordinary updates\.)",
    CBLOCK, "tensor-update implementation")
m = m.replace(
    " * Only the aligned coefficient m11 and cross coefficient m12 act\n"
    " * on this component. m22 is retained for tensor storage/validation.",
    " * Diagonal entries remain m11. Cross terms use reciprocal edge pairs\n"
    " * and normalized correlations derived from m11, m12, and m22.")
m, count = re.subn(r"update_tensor_edge\(g,\s*t,", "update_tensor_edge(g, mat, t,", m)
if count != 2:
    raise SystemExit(f"Expected two tensor-update call sites; found {count}. No files changed.")

# Analytic uniform-medium oracle needs all reciprocal neighbors present.
start = t.index("static void test_analytic_update(")
end = t.index("static void test_diagonal_invariance(", start)
part = t[start:end]
part = once(part,
    r"    REQUIRE\(materials_set_component_tensor\(\s*s\.m, s\.g, FIELD_EX, i, j, 0,\s*mxx, mxy, myy\) == 0\);\s*REQUIRE\(materials_set_component_tensor\(\s*s\.m, s\.g, FIELD_EY, i, j, 0,\s*myy, mxy, mxx\) == 0\);",
    """    /* Uniform patch supplies all four reciprocal neighbors. */
    for (size_t ii = i - 1; ii <= i + 1; ii++) {
        for (size_t jj = j - 1; jj <= j + 1; jj++) {
            REQUIRE(materials_set_component_tensor(
                s.m, s.g, FIELD_EX, ii, jj, 0,
                mxx, mxy, myy) == 0);
            REQUIRE(materials_set_component_tensor(
                s.m, s.g, FIELD_EY, ii, jj, 0,
                myy, mxy, mxx) == 0);
        }
    }""", "analytic tensor fixture")
t = t[:start] + part + t[end:]

# Add a one-cell tensor halo; checked Fourier-oracle points are unchanged.
start = t.index("static void test_periodic_tensor_seams(")
end = t.index("static void test_periodic_tensor_seam_suite(", start)
part = t[start:end]
pattern = r"    for \(size_t i = ibegin; i < iend; i\+\+\)\s*\{\s*for \(size_t j = jbegin; j < jend; j\+\+\)\s*\{"
replacement = """    size_t ti0 = (periodic_axes & 1u) ? ibegin : ibegin - 1;
    size_t ti1 = (periodic_axes & 1u) ? iend : iend + 1;
    size_t tj0 = (periodic_axes & 2u) ? jbegin : jbegin - 1;
    size_t tj1 = (periodic_axes & 2u) ? jend : jend + 1;

    for (size_t i = ti0; i < ti1; i++) {
        for (size_t j = tj0; j < tj1; j++) {"""
part, count = re.subn(pattern, replacement, part, count=1)
if count != 1:
    raise SystemExit("Could not find seam assignment loop. No files changed.")
t = t[:start] + part + t[end:]

# Symmetry is now an acceptance condition for every variable pattern.
start = t.index("static void test_variable_tensor_operator(")
end = t.index("static void test_variable_tensor_evolution(", start)
part = t[start:end]
part = once(part,
    r"    /\*\s*\* Uniform coefficients and diagonal-only variable coefficients.*?if \(pattern == 0u \|\| coupling_sign == 0\.0\)\s*\{\s*CHECK\(max_asymmetry <=\s*1e-11 \* fmax\(1\.0, max_entry\)\);\s*\}",
    """    /* Reciprocal assembly must be symmetric for every pattern. */
    CHECK(max_asymmetry <= 1e-11 * fmax(1.0, max_entry));""",
    "variable-operator symmetry assertion")
t = t[:start] + part + t[end:]

# Define the changed isolated-edge behavior explicitly.
extra = r"""
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

"""
t = once(t, r"int main\(void\)\s*\{", extra + "int main(void) {", "main function")
t = once(t, r"    test_variable_tensor_suite\(\);", "    test_variable_tensor_suite();\n    test_isolated_tensor_edge();", "variable-suite call")

stamp = time.strftime("%Y%m%d-%H%M%S")
for path in (mp, tp):
    backup = path.with_name(path.name + ".before-pair-fix-" + stamp)
    if backup.exists():
        raise SystemExit("Backup already exists; no files changed.")
for path in (mp, tp):
    shutil.copy2(path, path.with_name(path.name + ".before-pair-fix-" + stamp))
mp.write_text(m)
tp.write_text(t)
print("Patched materials.c and materials_tensor_unit.c; timestamped backups retained.")
print("No headers or CMake configuration changed; no Git commands were run.")
print("This candidate still requires both focused builds and regression testing.")
