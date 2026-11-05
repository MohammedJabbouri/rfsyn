import argparse, hashlib, re, shutil
from pathlib import Path
from datetime import datetime, timezone

ASSIGN = "int materials_assign_box(materials_t *mat, const fdtd_grid_t *g,\n    size_t i0, size_t i1, size_t j0, size_t j1,\n    size_t k0, size_t k1, material_id_t id)\n{\n    if (!box_valid(mat, g, i0, i1, j0, j1, k0, k1) ||\n        (size_t)id >= mat->table_count) return -1;\n\n    /* Half-open boxes are geometric ranges, not wrapped index lists.\n       Periodic endpoint storage must never alias plane zero here. */\n    if ((g->periodic_axes & 1u) && i1 > g->nx) i1 = g->nx;\n    if ((g->periodic_axes & 2u) && j1 > g->ny) j1 = g->ny;\n    if ((g->periodic_axes & 4u) && k1 > g->nz) k1 = g->nz;\n    if (i0 >= i1 || j0 >= j1 || k0 >= k1) return -1;\n\n    for (size_t i = i0; i < i1; i++) {\n        for (size_t j = j0; j < j1; j++) {\n            size_t idx = (i * (g->ny + 1) + j) * (g->nz + 1) + k0;\n            memset(mat->ids + idx, id, k1 - k0);\n        }\n    }\n    size_t kept = 0;\n    for (size_t n = 0; n < mat->override_count; n++) {\n        size_t idx = mat->overrides[n].key / 3;\n        size_t k = idx % (g->nz + 1);\n        size_t row = idx / (g->nz + 1);\n        size_t j = row % (g->ny + 1);\n        size_t i = row / (g->ny + 1);\n        if (i >= i0 && i < i1 && j >= j0 && j < j1 &&\n            k >= k0 && k < k1) continue;\n        mat->overrides[kept++] = mat->overrides[n];\n    }\n    mat->override_count = kept;\n    return 0;\n}\n"
GEOMETRY = "static inline int xt_medium(xt_sim *s, double eps_r)\n{\n    if (!s || !s->g || !s->m) return -1;\n    const fdtd_grid_t *g = s->g;\n    size_t ix = g->nx + ((g->periodic_axes & 1u) ? 0u : 1u);\n    size_t jy = g->ny + ((g->periodic_axes & 2u) ? 0u : 1u);\n    size_t kz = g->nz + ((g->periodic_axes & 4u) ? 0u : 1u);\n    return materials_set_box(s->m, g, 0, ix, 0, jy, 0, kz, eps_r, 0.0);\n}\n\nstatic inline int xt_interface(xt_sim *s, size_t i_interface,\n    double eps_left, double eps_right)\n{\n    if (!s || !s->g || !s->m || i_interface == 0 ||\n        i_interface >= s->g->nx) return -1;\n    if (xt_medium(s, eps_left) != 0) return -1;\n    const fdtd_grid_t *g = s->g;\n    size_t ix = g->nx + ((g->periodic_axes & 1u) ? 0u : 1u);\n    size_t jy = g->ny + ((g->periodic_axes & 2u) ? 0u : 1u);\n    size_t kz = g->nz + ((g->periodic_axes & 4u) ? 0u : 1u);\n    if (materials_set_box(s->m, g, i_interface, ix, 0, jy, 0, kz,\n                         eps_right, 0.0) != 0) return -1;\n\n    double eps_tangent = 0.5 * (eps_left + eps_right);\n    size_t planes[2] = {i_interface, 0};\n    size_t count = (g->periodic_axes & 1u) ? 2u : 1u;\n    for (size_t p = 0; p < count; p++) {\n        for (size_t j = 0; j < g->ny; j++) {\n            for (size_t k = 0; k < g->nz; k++) {\n                if (materials_set_component(s->m, g, FIELD_EY,\n                    planes[p], j, k, eps_tangent, 0.0) != 0) return -1;\n                if (materials_set_component(s->m, g, FIELD_EZ,\n                    planes[p], j, k, eps_tangent, 0.0) != 0) return -1;\n            }\n        }\n    }\n    return 0;\n}\n"
TEST = '#include "../engine/fdtd/extended_test.h"\n\nstatic int errors = 0;\n#define CHECK(x) do { if (!(x)) { fprintf(stderr,"FAIL line %d: %s\\n",__LINE__,#x); errors++; } } while (0)\n\nstatic void expect_epsilon(xt_sim *s, field_component_t component,\n    size_t i, size_t j, size_t k, double expected)\n{\n    material_coeff_t c;\n    int status = materials_get_coefficients(s->m,s->g,component,i,j,k,&c);\n    CHECK(status == 0);\n    if (status != 0) return;\n    CHECK(isfinite(c.Cb) && c.Cb > 0.0);\n    CHECK(fabs(c.Ca-1.0) < 1e-14);\n    double measured = s->g->dt/(s->g->dx*c.Cb*EPS0);\n    if (fabs(measured-expected) > 1e-11*fmax(1.0,expected)) {\n        fprintf(stderr,"coefficient mismatch component=%d (%zu,%zu,%zu): %.17g expected %.17g\\n",\n                (int)component,i,j,k,measured,expected);\n        errors++;\n    }\n}\n\nstatic void interface_map(unsigned axes)\n{\n    xt_sim s = {0};\n    int status=xt_create(&s,8,6,4,1e-3,axes,0);\n    CHECK(status==0); if(status!=0)return;\n    status=xt_interface(&s,4,1.0,25.0);\n    CHECK(status==0);\n    if (status==0) {\n        for(size_t i=0;i<8;i++) for(size_t j=0;j<6;j++) for(size_t k=0;k<4;k++) {\n            double base=i<4?1.0:25.0;\n            double tangent=(i==4 || ((axes&1u)&&i==0))?13.0:base;\n            expect_epsilon(&s,FIELD_EX,i,j,k,base);\n            expect_epsilon(&s,FIELD_EY,i,j,k,tangent);\n            expect_epsilon(&s,FIELD_EZ,i,j,k,tangent);\n        }\n        if (!(axes&1u)) {\n            expect_epsilon(&s,FIELD_EY,8,0,0,25.0);\n            expect_epsilon(&s,FIELD_EZ,8,0,0,25.0);\n        }\n    }\n    xt_destroy(&s);\n}\n\nstatic void assignment_boxes(void)\n{\n    xt_sim s={0};\n    int status=xt_create(&s,8,6,4,1e-3,7u,0);\n    CHECK(status==0);if(status!=0)return;\n    material_id_t dielectric;\n    CHECK(materials_register(s.m,25.0,0.0,&dielectric)==0);\n    CHECK(materials_set_component(s.m,s.g,FIELD_EY,0,0,0,7.0,0.0)==0);\n    CHECK(materials_set_component(s.m,s.g,FIELD_EY,4,3,2,9.0,0.0)==0);\n    CHECK(materials_assign_box(s.m,s.g,4,9,3,7,2,5,dielectric)==0);\n    for(size_t i=0;i<8;i++) for(size_t j=0;j<6;j++) for(size_t k=0;k<4;k++) {\n        double base=(i>=4&&j>=3&&k>=2)?25.0:1.0;\n        expect_epsilon(&s,FIELD_EX,i,j,k,base);\n        expect_epsilon(&s,FIELD_EY,i,j,k,(i==0&&j==0&&k==0)?7.0:base);\n        expect_epsilon(&s,FIELD_EZ,i,j,k,base);\n    }\n    CHECK(materials_assign_box(s.m,s.g,8,9,0,1,0,1,dielectric)!=0);\n    expect_epsilon(&s,FIELD_EY,0,0,0,7.0);\n    CHECK(materials_assign_box(s.m,s.g,0,9,0,7,0,5,dielectric)==0);\n    for(size_t i=0;i<8;i++) for(size_t j=0;j<6;j++) for(size_t k=0;k<4;k++)\n        for(int c=0;c<3;c++) expect_epsilon(&s,(field_component_t)c,i,j,k,25.0);\n    xt_destroy(&s);\n}\n\nint main(void)\n{\n    interface_map(7u);\n    interface_map(6u);\n    interface_map(0u);\n    assignment_boxes();\n    if(errors)fprintf(stderr,"%d periodic material checks failed\\n",errors);\n    else printf("PASS: periodic box assignment, overrides, interface and seam coefficients; nonperiodic control\\n");\n    return errors?1:0;\n}\n'
CMAKE = r"""

if(BUILD_TESTING)
  if(NOT TARGET periodic_materials_unit)
    add_executable(periodic_materials_unit tests/periodic_materials_unit.c)
    target_link_libraries(periodic_materials_unit PRIVATE rfsyn_fdtd ${RFSYN_MATH_LIBRARIES})
    rfsyn_apply_build_options(periodic_materials_unit)
    add_test(NAME periodic_materials_unit COMMAND $<TARGET_FILE:periodic_materials_unit>)
    set_tests_properties(periodic_materials_unit PROPERTIES LABELS "unit;fdtd;periodic")
  endif()
endif()
"""


def blob(text):
    data = text.encode()
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()


def main():
    parser = argparse.ArgumentParser(
        description="Prepare/apply the narrowly scoped periodic material fix. Never commits."
    )
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Write files after validating all inputs; default is dry run",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    material = root / "engine/fdtd/materials.c"
    header = root / "engine/fdtd/extended_test.h"
    cmake = root / "CMakeLists.txt"
    test = root / "tests/periodic_materials_unit.c"
    original = {p: p.read_text() for p in [material, header, cmake]}
    expected = {
        material: "e3ad56e15e6c17510bba0b575b1e539370da0375",
        header: "2ea0743475fce6aaba7e377c5a6c18b8d82da0f7",
    }
    for p, sha in expected.items():
        if blob(original[p]) != sha:
            raise SystemExit(
                f"STOP: {p} differs from reviewed source. No files changed. Review the local file before patching."
            )
    if test.exists():
        raise SystemExit(f"STOP: {test} already exists; no files changed")
    if "periodic_materials_unit" in original[cmake]:
        raise SystemExit(
            "STOP: existing periodic_materials_unit CMake entry; no files changed"
        )
    a = original[material]
    start = a.index("int materials_assign_box(")
    end = a.index("int materials_set_box(", start)
    material_new = a[:start] + ASSIGN + "\n" + a[end:]
    h = original[header]
    start = h.index("static inline int xt_medium(")
    end = h.index("static inline double *xt_e(", start)
    header_new = h[:start] + GEOMETRY + "\n" + h[end:]
    changes = {
        material: material_new,
        header: header_new,
        cmake: original[cmake] + "\n" + CMAKE,
        test: TEST,
    }
    import difflib

    for p, new in changes.items():
        old = original.get(p, "")
        print(
            "".join(
                difflib.unified_diff(
                    old.splitlines(True),
                    new.splitlines(True),
                    fromfile=str(p.relative_to(root)),
                    tofile=str(p.relative_to(root)),
                )
            ),
            end="",
        )
    if not args.apply:
        print(
            "\nDRY RUN: no files changed. Re-run with --apply after reviewing the diff."
        )
        return
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    for p in original:
        shutil.copy2(p, p.with_name(p.name + ".before-periodic-" + timestamp))
    for p, new in changes.items():
        p.write_text(new)
    print(
        "\nApplied. Existing files have timestamped backups. Build and run tests; nothing was committed."
    )


if __name__ == "__main__":
    main()
