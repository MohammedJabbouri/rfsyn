cmake_minimum_required(VERSION 3.16)

project(rfsyn LANGUAGES C)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(
  RFSYN_SANITIZE
  "Enable address and undefined-behavior sanitizers"
  OFF
)

option(
  RFSYN_ENABLE_OPENMP
  "Enable OpenMP when available"
  ON
)

option(
  RFSYN_BUILD_EXTENDED_VALIDATORS
  "Build extended validators and tilted-interface measurements"
  ON
)

set(
  RFSYN_VALIDATION_REFINEMENT
  "1"
  CACHE STRING
  "Supported validator refinement: 1, 2, or 4"
)

set_property(
  CACHE RFSYN_VALIDATION_REFINEMENT
  PROPERTY STRINGS 1 2 4
)

if(NOT "${RFSYN_VALIDATION_REFINEMENT}" MATCHES "^(1|2|4)$")
  message(FATAL_ERROR "RFSYN_VALIDATION_REFINEMENT must be 1, 2, or 4")
endif()

set(
  RFSYN_STABILITY_STEPS
  "100000"
  CACHE STRING
  "Number of steps per stability case"
)

if(NOT "${RFSYN_STABILITY_STEPS}" MATCHES "^[1-9][0-9]*$")
  message(FATAL_ERROR "RFSYN_STABILITY_STEPS must be a positive integer")
endif()

if(RFSYN_STABILITY_STEPS LESS 10000 OR
   RFSYN_STABILITY_STEPS GREATER 10000000)
  message(
    FATAL_ERROR
    "RFSYN_STABILITY_STEPS must be between 10000 and 10000000"
  )
endif()

set(
  RFSYN_VALIDATION_TIMEOUT
  "900"
  CACHE STRING
  "Timeout in seconds per baseline validator"
)

set(
  RFSYN_EXTENDED_TIMEOUT
  "7200"
  CACHE STRING
  "Timeout in seconds per extended validator or measurement"
)

foreach(timeout_var IN ITEMS
    RFSYN_VALIDATION_TIMEOUT
    RFSYN_EXTENDED_TIMEOUT)
  if(NOT "${${timeout_var}}" MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "${timeout_var} must be a positive integer")
  endif()
endforeach()

# ---------------------------------------------------------------------------
# Shared build options
# ---------------------------------------------------------------------------

add_library(rfsyn_build_options INTERFACE)

if(RFSYN_SANITIZE)
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang" AND NOT MSVC)
    target_compile_options(
      rfsyn_build_options
      INTERFACE
        -fsanitize=address,undefined
        -fno-omit-frame-pointer
        -g
    )

    target_link_options(
      rfsyn_build_options
      INTERFACE
        -fsanitize=address,undefined
    )
  else()
    message(
      FATAL_ERROR
      "Sanitizers require a supported GCC or Clang toolchain"
    )
  endif()
endif()

set(RFSYN_MATH_LIBRARIES)

if(NOT MSVC)
  list(APPEND RFSYN_MATH_LIBRARIES m)
endif()

function(rfsyn_apply_build_options target_name)
  target_link_libraries(
    ${target_name}
    PRIVATE
      rfsyn_build_options
  )
endfunction()

# ---------------------------------------------------------------------------
# FDTD library
# ---------------------------------------------------------------------------

add_library(
  rfsyn_fdtd STATIC
  engine/fdtd/grid.c
  engine/fdtd/cpml.c
  engine/fdtd/materials.c
  engine/fdtd/source.c
  engine/fdtd/probe.c
  engine/fdtd/simulate.c
  engine/fdtd/subpixel_te.c
)

if(RFSYN_MATH_LIBRARIES)
  target_link_libraries(
    rfsyn_fdtd
    PUBLIC
      ${RFSYN_MATH_LIBRARIES}
  )
endif()

rfsyn_apply_build_options(rfsyn_fdtd)

if(RFSYN_ENABLE_OPENMP)
  find_package(OpenMP QUIET COMPONENTS C)

  if(OpenMP_C_FOUND)
    target_link_libraries(
      rfsyn_fdtd
      PUBLIC
        OpenMP::OpenMP_C
    )

    target_compile_definitions(
      rfsyn_fdtd
      PUBLIC
        RFSYN_HAS_OPENMP=1
    )

    message(STATUS "rfsyn: OpenMP enabled")
  else()
    message(STATUS "rfsyn: OpenMP unavailable; FDTD is single-threaded")
  endif()
else()
  message(STATUS "rfsyn: OpenMP disabled")
endif()

# ---------------------------------------------------------------------------
# Core library
# ---------------------------------------------------------------------------

add_library(
  rfsyn_core STATIC
  engine/core/chain.c
  engine/core/chains.c
  engine/core/rf_signal.c
  engine/core/rng.c
)

if(RFSYN_MATH_LIBRARIES)
  target_link_libraries(
    rfsyn_core
    PUBLIC
      ${RFSYN_MATH_LIBRARIES}
  )
endif()

rfsyn_apply_build_options(rfsyn_core)

# ---------------------------------------------------------------------------
# JSON library
# ---------------------------------------------------------------------------

add_library(
  rfsyn_cjson STATIC
  cli/cJSON.c
)

target_include_directories(
  rfsyn_cjson
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/cli"
)

rfsyn_apply_build_options(rfsyn_cjson)

# ---------------------------------------------------------------------------
# World library
# ---------------------------------------------------------------------------

add_library(
  rfsyn_world STATIC
  engine/world/world.c
  engine/world/material_catalog.c
  cli/material_config.c
)

target_include_directories(
  rfsyn_world
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(
  rfsyn_world
  PUBLIC
    rfsyn_cjson
    ${RFSYN_MATH_LIBRARIES}
)

rfsyn_apply_build_options(rfsyn_world)

# ---------------------------------------------------------------------------
# Noise library
# ---------------------------------------------------------------------------

add_library(
  rfsyn_noise STATIC
  engine/noise/awgn.c
)

target_link_libraries(
  rfsyn_noise
  PUBLIC
    rfsyn_core
    ${RFSYN_MATH_LIBRARIES}
)

rfsyn_apply_build_options(rfsyn_noise)

# ---------------------------------------------------------------------------
# Application
# ---------------------------------------------------------------------------

add_executable(
  rfsyn
  cli/main.c
  cli/settings.c
  cli/config.c
  cli/presets.c
  cli/fdtd_stage.c
)

target_link_libraries(
  rfsyn
  PRIVATE
    rfsyn_fdtd
    rfsyn_core
    rfsyn_noise
    rfsyn_world
    ${RFSYN_MATH_LIBRARIES}
)

rfsyn_apply_build_options(rfsyn)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

include(CTest)

if(BUILD_TESTING)



  # -------------------------------------------------------------------------
  # Unit-test registration
  # -------------------------------------------------------------------------

  function(rfsyn_add_unit_test name)
    add_executable(
      ${name}
      "tests/${name}.c"
    )

    target_link_libraries(
      ${name}
      PRIVATE
        ${ARGN}
        ${RFSYN_MATH_LIBRARIES}
    )

    rfsyn_apply_build_options(${name})

    add_test(
      NAME ${name}
      COMMAND $<TARGET_FILE:${name}>
    )

    set_tests_properties(
      ${name}
      PROPERTIES
        LABELS "unit"
    )
  endfunction()

  rfsyn_add_unit_test(world_unit rfsyn_world)
  rfsyn_add_unit_test(material_catalog_unit rfsyn_world)
  rfsyn_add_unit_test(material_config_unit rfsyn_world)
  rfsyn_add_unit_test(simulate_unit rfsyn_fdtd)
  rfsyn_add_unit_test(materials_unit rfsyn_fdtd)
  rfsyn_add_unit_test(rf_signal_unit rfsyn_core)
  rfsyn_add_unit_test(subpixel_te_unit rfsyn_fdtd)

  # -------------------------------------------------------------------------
  # CLI smoke test
  # -------------------------------------------------------------------------

  add_executable(
    cli_smoke
    tests/cli_smoke.c
  )

  target_link_libraries(
    cli_smoke
    PRIVATE
      rfsyn_core
      rfsyn_noise
      rfsyn_fdtd
      rfsyn_world
      ${RFSYN_MATH_LIBRARIES}
  )

  rfsyn_apply_build_options(cli_smoke)

  add_dependencies(cli_smoke rfsyn)

  add_test(
    NAME cli_smoke
    COMMAND
      $<TARGET_FILE:cli_smoke>
      $<TARGET_FILE:rfsyn>
  )

  set_tests_properties(
    cli_smoke
    PROPERTIES
      LABELS "smoke;cli"
  )

  # -------------------------------------------------------------------------
  # Validator and measurement registration
  #
  # Additional arguments are forwarded to the executable.
  # -------------------------------------------------------------------------

  function(rfsyn_add_validator name group timeout_seconds)
    set(
      validator_source
      "${CMAKE_CURRENT_SOURCE_DIR}/engine/fdtd/${name}.c"
    )

    if(NOT EXISTS "${validator_source}")
      message(
        FATAL_ERROR
        "Missing validator source: ${validator_source}"
      )
    endif()

    add_executable(
      ${name}
      "${validator_source}"
    )

    target_link_libraries(
      ${name}
      PRIVATE
        rfsyn_fdtd
        ${RFSYN_MATH_LIBRARIES}
    )

    target_compile_definitions(
      ${name}
      PRIVATE
        RFSYN_VALIDATION_REFINEMENT=${RFSYN_VALIDATION_REFINEMENT}
    )

    rfsyn_apply_build_options(${name})

    set(
      validator_workdir
      "${CMAKE_CURRENT_BINARY_DIR}/validation/${name}"
    )

    file(MAKE_DIRECTORY "${validator_workdir}")

    add_test(
      NAME ${name}
      COMMAND $<TARGET_FILE:${name}> ${ARGN}
    )

    if("${group}" STREQUAL "measurement")
      set(validator_labels "fdtd;measurement;tilted")
    else()
      set(validator_labels "fdtd;validation;${group}")
    endif()

    set_tests_properties(
      ${name}
      PROPERTIES
        WORKING_DIRECTORY "${validator_workdir}"
        LABELS "${validator_labels}"
        TIMEOUT "${timeout_seconds}"
    )
  endfunction()

  # -------------------------------------------------------------------------
  # Existing six baseline validators
  # -------------------------------------------------------------------------

  set(
    RFSYN_FDTD_VALIDATORS
    validate_cpml
    validate_cpml_margin
    validate_cpml_planewave
    validate_materials_fresnel
    validate_probe_dft
    validate_materials_energy
  )

  foreach(name IN LISTS RFSYN_FDTD_VALIDATORS)
    rfsyn_add_validator(
      ${name}
      baseline
      "${RFSYN_VALIDATION_TIMEOUT}"
    )
  endforeach()

  # -------------------------------------------------------------------------
  # Extended validators and tilted-interface measurement
  # -------------------------------------------------------------------------

  if(RFSYN_BUILD_EXTENDED_VALIDATORS)

    set(
      RFSYN_REQUIRED_EXTENDED_FILES
      engine/fdtd/extended_test.h
      engine/fdtd/validate_cpml_oblique.c
      engine/fdtd/validate_materials_high_contrast.c
      engine/fdtd/validate_long_term_stability.c
      engine/fdtd/validate_materials_tilted.c
    )

    set(RFSYN_MISSING_EXTENDED_FILES)

    foreach(relative_path IN LISTS RFSYN_REQUIRED_EXTENDED_FILES)
      if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${relative_path}")
        list(
          APPEND
          RFSYN_MISSING_EXTENDED_FILES
          "${relative_path}"
        )
      endif()
    endforeach()

    if(RFSYN_MISSING_EXTENDED_FILES)
      string(
        REPLACE ";" "\n  "
        missing_files_text
        "${RFSYN_MISSING_EXTENDED_FILES}"
      )

      message(
        FATAL_ERROR
        "Missing extended-validator files:\n"
        "  ${missing_files_text}\n"
        "Add these files or configure with "
        "-DRFSYN_BUILD_EXTENDED_VALIDATORS=OFF."
      )
    endif()

    # Declaration checks only; these do not verify implementation behavior.
    file(
      READ
      "${CMAKE_CURRENT_SOURCE_DIR}/engine/fdtd/grid.h"
      grid_header
    )

    if(NOT grid_header MATCHES "fdtd_grid_set_periodic")
      message(
        FATAL_ERROR
        "Extended validators require the periodic-grid accuracy rewrite "
        "and fdtd_grid_set_periodic()."
      )
    endif()

    file(
      READ
      "${CMAKE_CURRENT_SOURCE_DIR}/engine/fdtd/materials.h"
      materials_header
    )

    foreach(required_api IN ITEMS
        materials_register
        materials_assign_box
        materials_set_component
        materials_get_coefficients)
      if(NOT materials_header MATCHES "${required_api}")
        message(
          FATAL_ERROR
          "Required materials API missing: ${required_api}"
        )
      endif()
    endforeach()

    rfsyn_add_validator(
      validate_cpml_oblique
      extended
      "${RFSYN_EXTENDED_TIMEOUT}"
    )

    rfsyn_add_validator(
      validate_materials_high_contrast
      extended
      "${RFSYN_EXTENDED_TIMEOUT}"
      --full
    )

    rfsyn_add_validator(
      validate_long_term_stability
      extended
      "${RFSYN_EXTENDED_TIMEOUT}"
      --steps
      "${RFSYN_STABILITY_STEPS}"
    )

    rfsyn_add_validator(
      validate_materials_tilted
      measurement
      "${RFSYN_EXTENDED_TIMEOUT}"
    )

  endif()

elseif(RFSYN_BUILD_EXTENDED_VALIDATORS)
  message(
    STATUS
    "rfsyn: BUILD_TESTING=OFF; validators will not be built"
  )
endif()

# ---------------------------------------------------------------------------
# Configuration summary
# ---------------------------------------------------------------------------

message(STATUS "rfsyn: C compiler = ${CMAKE_C_COMPILER}")
message(STATUS "rfsyn: tests = ${BUILD_TESTING}")
message(STATUS "rfsyn: refinement = ${RFSYN_VALIDATION_REFINEMENT}")
message(STATUS "rfsyn: extended validators = ${RFSYN_BUILD_EXTENDED_VALIDATORS}")
message(STATUS "rfsyn: stability steps per case = ${RFSYN_STABILITY_STEPS}")
message(STATUS "rfsyn: sanitizers = ${RFSYN_SANITIZE}")