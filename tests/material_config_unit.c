#include <stdio.h>

#include <math.h>

#include <string.h>

#include "../cli/material_config.h"

#include "cJSON.h"

#define RC_OK 0
#define RC_PARSE - 2

#define CHECK(condition) do { if (!(condition)) { fprintf(stderr, "FAIL at line %d: %s\n", __LINE__, #condition); goto fail; } } while (0)

static material_catalog_t * g_cat = NULL;
static cJSON * g_root = NULL;
static int g_failed = 0;

static void reset_scene(const char * json_text) {
  material_catalog_destroy(g_cat);
  cJSON_Delete(g_root);
  g_cat = material_catalog_create_default();
  g_root = json_text ? cJSON_Parse(json_text) : NULL;
}

int main(void) {
  g_failed = 1;

  register int dummy = 0;
  (void) dummy;

  reset_scene(NULL);
  {
    material_config_result_t r;
    CHECK(material_config_load(NULL, g_cat, & r, stderr) == RC_OK);
    CHECK(r.enabled == 0 && r.allow_proxy == 0);
    CHECK(r.evaluation_frequency_hz == 1e9);
    CHECK(r.materials_loaded == 0 && r.colors_seen == 0);
  }

  reset_scene("{}");
  {
    material_config_result_t r;
    CHECK(material_config_load(g_root, g_cat, & r, stderr) == RC_OK);
    CHECK(r.materials_loaded == 0 && r.colors_seen == 0);
  }

  reset_scene(
    "{\"materials\":{\"my_wall\":{\"model\":\"constant\",\"eps_r\":4.0,"
    "\"sigma_s_per_m\":0.01,\"reference_frequency_hz\":1000000000,"
    "\"description\":\"Test-only definition\",\"source\":\"Unit test; not measured\"}}}");
  {
    material_config_result_t r;
    catalog_properties_t p;
    CHECK(material_config_load(g_root, g_cat, & r, stderr) == RC_OK);
    CHECK(r.materials_loaded == 1);
    CHECK(material_catalog_evaluate(g_cat, "my_wall", 1e9, 0, & p) == CATALOG_OK);
    CHECK(p.eps_r == 4.0 && p.sigma_s_per_m == 0.01);
    CHECK(p.status == CATALOG_USER_DEFINED);
  }

  reset_scene(
    "{\"materials\":{\"concrete_normal\":{\"model\":\"constant\",\"eps_r\":4.0,"
    "\"sigma_s_per_m\":0.01}}}");
  CHECK(material_config_load(g_root, g_cat, & (material_config_result_t) {
    0
  }, stderr) != RC_OK);

  reset_scene(
    "{\"materials\":{\"concrete_normal\":{\"model\":\"constant\",\"eps_r\":4.0,"
    "\"sigma_s_per_m\":0.01,\"override\":true,\"source\":\"Unit-test override\"}}}");
  {
    material_config_result_t r;
    catalog_material_id_t id, ref;
    catalog_properties_t p;
    CHECK(material_catalog_find(g_cat, "concrete_normal", & ref) == CATALOG_OK);
    CHECK(material_config_load(g_root, g_cat, & r, stderr) == RC_OK);
    CHECK(material_catalog_find(g_cat, "concrete_normal", & id) == CATALOG_OK && id == ref);
    CHECK(material_catalog_evaluate(g_cat, "concrete_normal", 1e9, 0, & p) == CATALOG_OK);
    CHECK(p.eps_r == 4.0 && p.status == CATALOG_USER_DEFINED);
  }

  reset_scene(
    "{\"materials\":{\"bad\":{\"eps_r\":4.0,\"sigma_s_per_m\":0.01}}}");
  CHECK(material_config_load(g_root, g_cat, & (material_config_result_t) {
    0
  }, stderr) != RC_OK);

  reset_scene(
    "{\"materials\":{\"bad\":{\"model\":\"constant\",\"eps_r\":0.5,"
    "\"sigma_s_per_m\":0.01}}}");
  CHECK(material_config_load(g_root, g_cat, & (material_config_result_t) {
    0
  }, stderr) != RC_OK);

  reset_scene(
    "{\"materials\":{\"air\":{\"model\":\"constant\",\"eps_r\":1.0,"
    "\"sigma_s_per_m\":0.0,\"override\":true}}}");
  CHECK(material_config_load(g_root, g_cat, & (material_config_result_t) {
    0
  }, stderr) != RC_OK);

  reset_scene(
    "{\"material_catalog\":{\"enabled\":false,\"evaluation_frequency_hz\":-5}}");
  CHECK(material_config_load(g_root, g_cat, & (material_config_result_t) {
    0
  }, stderr) != RC_OK);

  reset_scene(
    "{\"material_catalog\":{\"enabled\":true,\"evaluation_frequency_hz\":1000000000,"
    "\"allow_proxy\":true},"
    "\"colors\":{\"#808080\":\"concrete\",\"#F0E8D8\":\"plasterboard\"}}");
  {
    material_config_result_t r;
    char name[MATERIAL_CATALOG_NAME_SIZE];
    CHECK(material_config_load(g_root, g_cat, & r, stderr) == RC_OK);
    CHECK(r.enabled == 1 && r.allow_proxy == 1);
    CHECK(r.colors_seen == 2);
    CHECK(material_catalog_lookup_color(g_cat, (color_rgb_t) {
        0x80,
        0x80,
        0x80
      },
      name, sizeof(name)) == CATALOG_OK);
    CHECK(strcmp(name, "concrete_normal") == 0);
    CHECK(material_catalog_lookup_color(g_cat, (color_rgb_t) {
        0xF0,
        0xE8,
        0xD8
      },
      name, sizeof(name)) == CATALOG_OK);
    CHECK(strcmp(name, "plasterboard") == 0);
  }

  reset_scene(
    "{\"material_catalog\":{\"enabled\":true},"
    "\"colors\":{\"#808080\":\"material_that_does_not_exist\"}}");
  CHECK(material_config_load(g_root, g_cat, & (material_config_result_t) {
    0
  }, stderr) != RC_OK);

  reset_scene(
    "{\"material_catalog\":{\"enabled\":true},"
    "\"colors\":{\"red\":\"concrete_normal\"}}");
  CHECK(material_config_load(g_root, g_cat, & (material_config_result_t) {
    0
  }, stderr) != RC_OK);

  reset_scene(
    "{\"material_catalog\":{\"enabled\":false},"
    "\"colors\":{\"red\":\"concrete_normal\"}}");
  {
    material_config_result_t r;
    CHECK(material_config_load(g_root, g_cat, & r, stderr) == RC_OK);
    CHECK(r.enabled == 0 && r.colors_seen == 0);
  }

  {
    cJSON * bad = cJSON_Parse("[1,2,3]");
    CHECK(bad != NULL);
    CHECK(material_config_load(bad, g_cat, & (material_config_result_t) {
      0
    }, stderr) == RC_PARSE);
    cJSON_Delete(bad);
  }

  /* Case 15: materials must be an object. */
  reset_scene("{\"materials\":42}");
  CHECK(material_config_load(g_root, g_cat, & (material_config_result_t) {
    0
  }, stderr) != RC_OK);

  /* Case 16: alias names cannot be defined directly. */
  reset_scene(
    "{\"materials\":{\"vacuum\":{\"model\":\"constant\",\"eps_r\":1.0,"
    "\"sigma_s_per_m\":0.0}}}");
  CHECK(material_config_load(g_root, g_cat, & (material_config_result_t) {
    0
  }, stderr) != RC_OK);

  reset_scene(
    "{\"materials\":{\"asphalt\":{\"model\":\"constant\",\"eps_r\":4.0,"
    "\"sigma_s_per_m\":0.02,\"override\":true,\"source\":\"Illustration; not validated\"}}}");
  {
    catalog_properties_t p;
    CHECK(material_config_load(g_root, g_cat, & (material_config_result_t) {
      0
    }, stderr) == RC_OK);
    CHECK(material_catalog_evaluate(g_cat, "asphalt", 1e9, 0, & p) == CATALOG_OK);
    CHECK(p.eps_r == 4.0 && p.status == CATALOG_USER_DEFINED);
  }

  reset_scene("{\"material_catalog\":{\"evaluation_frequency_hz\":2400000000}}");
  {
    material_config_result_t r;
    CHECK(material_config_load(g_root, g_cat, & r, stderr) == RC_OK);
    CHECK(r.evaluation_frequency_hz == 2400000000.0);
  }

  printf("material_config_unit: PASS (32 cases)\n");
  {
    static
    const struct {
      const char * name;
      const char * json;
    }
    invalid_cases[] = {
      {
        "enabled must be boolean",
        "{\"material_catalog\":{\"enabled\":\"true\"}}"
      },
      {
        "allow_proxy must be boolean",
        "{\"material_catalog\":{\"allow_proxy\":1}}"
      },
      {
        "override must be boolean",
        "{\"materials\":{\"test_invalid\":{"
        "\"model\":\"constant\","
        "\"eps_r\":4.0,"
        "\"sigma_s_per_m\":0.01,"
        "\"override\":\"true\"}}}"
      },
      {
        "reference frequency must be numeric",
        "{\"materials\":{\"test_invalid\":{"
        "\"model\":\"constant\","
        "\"eps_r\":4.0,"
        "\"sigma_s_per_m\":0.01,"
        "\"reference_frequency_hz\":\"1000000000\"}}}"
      },
      {
        "reference frequency must be positive",
        "{\"materials\":{\"test_invalid\":{"
        "\"model\":\"constant\","
        "\"eps_r\":4.0,"
        "\"sigma_s_per_m\":0.01,"
        "\"reference_frequency_hz\":0}}}"
      },
      {
        "description must be string",
        "{\"materials\":{\"test_invalid\":{"
        "\"model\":\"constant\","
        "\"eps_r\":4.0,"
        "\"sigma_s_per_m\":0.01,"
        "\"description\":42}}}"
      },
      {
        "source must be string",
        "{\"materials\":{\"test_invalid\":{"
        "\"model\":\"constant\","
        "\"eps_r\":4.0,"
        "\"sigma_s_per_m\":0.01,"
        "\"source\":false}}}"
      },
      {
        "source must not be empty",
        "{\"materials\":{\"test_invalid\":{"
        "\"model\":\"constant\","
        "\"eps_r\":4.0,"
        "\"sigma_s_per_m\":0.01,"
        "\"source\":\"\"}}}"
      },
      {
        "evaluation frequency must be numeric",
        "{\"material_catalog\":{"
        "\"evaluation_frequency_hz\":\"1000000000\"}}"
      },
      {
        "unsupported schema version",
        "{\"material_catalog\":{\"schema_version\":2}}"
      },
      {
        "schema version must be numeric",
        "{\"material_catalog\":{\"schema_version\":\"1\"}}"
      },
      {
        "fractional schema version rejected",
        "{\"material_catalog\":{\"schema_version\":1.5}}"
      }
    };

    size_t case_count =
      sizeof(invalid_cases) / sizeof(invalid_cases[0]);

    for (size_t i = 0; i < case_count; i++) {
      material_config_result_t r;
      catalog_material_id_t invalid_id;
      catalog_properties_t air;
      int rc;

      reset_scene(invalid_cases[i].json);
      CHECK(g_cat != NULL);
      CHECK(g_root != NULL);

      rc = material_config_load(g_root, g_cat, & r, stderr);

      if (rc == RC_OK) {
        fprintf(stderr, "unexpected acceptance: %s\n",
          invalid_cases[i].name);
        goto fail;
      }

      CHECK(material_catalog_find(
          g_cat, "test_invalid", & invalid_id) ==
        CATALOG_NOT_FOUND);

      CHECK(material_catalog_evaluate(
        g_cat, "air", 1e9, 0, & air) == CATALOG_OK);
      CHECK(air.eps_r == 1.0);
      CHECK(air.sigma_s_per_m == 0.0);
    }
  }

  {
    material_config_result_t r;

    reset_scene(
      "{\"material_catalog\":{"
      "\"schema_version\":1,"
      "\"enabled\":true,"
      "\"allow_proxy\":true,"
      "\"evaluation_frequency_hz\":1000000000},"
      "\"materials\":{},"
      "\"colors\":{}}");

    CHECK(g_cat != NULL);
    CHECK(g_root != NULL);

    CHECK(material_config_load(
      g_root, g_cat, & r, stderr) == RC_OK);

    CHECK(r.enabled == 1);
    CHECK(r.allow_proxy == 1);
    CHECK(r.evaluation_frequency_hz == 1e9);
    CHECK(r.materials_loaded == 0);
    CHECK(r.colors_seen == 0);
  }

  {
    material_config_result_t r;
    catalog_material_id_t id;
    catalog_definition_t d;
    catalog_properties_t p;

    reset_scene(
      "{\"material_catalog\":{},"
      "\"materials\":{\"minimal_wall\":{"
      "\"model\":\"constant\","
      "\"eps_r\":4.0,"
      "\"sigma_s_per_m\":0.01}}}");

    CHECK(g_cat != NULL);
    CHECK(g_root != NULL);

    CHECK(material_config_load(
      g_root, g_cat, & r, stderr) == RC_OK);

    CHECK(r.enabled == 0);
    CHECK(r.allow_proxy == 0);
    CHECK(r.evaluation_frequency_hz == 1e9);
    CHECK(r.materials_loaded == 1);
    CHECK(r.colors_seen == 0);

    CHECK(material_catalog_find(
      g_cat, "minimal_wall", & id) == CATALOG_OK);

    CHECK(material_catalog_get(g_cat, id, & d) == CATALOG_OK);
    CHECK(d.reference_frequency_hz == 1e9);
    CHECK(strcmp(d.description, "") == 0);
    CHECK(strcmp(d.source, "User-supplied values") == 0);

    CHECK(material_catalog_evaluate(
      g_cat, "minimal_wall", 1e9, 0, & p) == CATALOG_OK);
    CHECK(p.eps_r == 4.0);
    CHECK(p.sigma_s_per_m == 0.01);
    CHECK(p.status == CATALOG_USER_DEFINED);
  }
  g_failed = 0;

  fail:
    material_catalog_destroy(g_cat);
  cJSON_Delete(g_root);
  return g_failed;
}