cmake_minimum_required(VERSION 3.16)

project(rfsyn LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(RFSYN_SANITIZE "Enable address and undefined-behavior sanitizers" OFF)
option(RFSYN_ENABLE_OPENMP "Enable OpenMP when available" ON)
option(RFSYN_BUILD_EXTENDED_VALIDATORS "Build oblique, high-contrast, and long-term validators" OFF)

set(
  RFSYN_VALIDATION_REFINEMENT
  "1"
  CACHE STRING
  "Normal-incidence validation refinement: 1, 2, or 4"
)

set_property(
  CACHE RFSYN_VALIDATION_REFINEMENT
  PROPERTY STRINGS 1 2 4
)

if(NOT "${RFSYN_VALIDATION_REFINEMENT}" MATCHES "^(1|2|4)$")
  message(FATAL_ERROR "RFSYN_VALIDATION_REFINEMENT must be 1, 2, or 4")
endif()

set(
  RFSYN_VALIDATION_TIMEOUT
  "900"
  CACHE STRING
  "Timeout in seconds for each standard FDTD validator"
)

set(
  RFSYN_EXTENDED_TIMEOUT
  "7200"
  CACHE STRING
  "Timeout in seconds for each extended FDTD validator"
)

foreach(timeout_var IN ITEMS RFSYN_VALIDATION_TIMEOUT RFSYN_EXTENDED_TIMEOUT)
  if(NOT "${${timeout_var}}" MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "${timeout_var} must be a positive integer")
  endif()
endforeach()

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
    message(FATAL_ERROR "RFSYN_SANITIZE requires a supported GCC or Clang toolchain")
  endif()
endif()

# MSVC provides these math functions through its runtime, not libm.
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

add_library(
  rfsyn_fdtd STATIC
  engine/fdtd/grid.c
  engine/fdtd/cpml.c
  engine/fdtd/materials.c
  engine/fdtd/source.c
  engine/fdtd/probe.c
  engine/fdtd/simulate.c
)

target_link_libraries(
  rfsyn_fdtd
  PUBLIC
    ${RFSYN_MATH_LIBRARIES}
)

rfsyn_apply_build_options(rfsyn_fdtd)

if(RFSYN_ENABLE_OPENMP)
  find_package(OpenMP QUIET COMPONENTS C)

  if(OpenMP_C_FOUND)
    target_link_libraries(rfsyn_fdtd PUBLIC OpenMP::OpenMP_C)
    target_compile_definitions(rfsyn_fdtd PUBLIC RFSYN_HAS_OPENMP=1)
    message(STATUS "rfsyn: OpenMP enabled")
  else()
    message(STATUS "rfsyn: OpenMP not found; using single-threaded FDTD")
  endif()
else()
  message(STATUS "rfsyn: OpenMP explicitly disabled")
endif()

add_library(
  rfsyn_core STATIC
  engine/core/chain.c
  engine/core/chains.c
  engine/core/rf_signal.c
  engine/core/rng.c
)

target_link_libraries(
  rfsyn_core
  PUBLIC
    ${RFSYN_MATH_LIBRARIES}
)

rfsyn_apply_build_options(rfsyn_core)

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

include(CTest)

if(BUILD_TESTING)

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

  # Building cli_smoke also builds the application it launches.
  add_dependencies(cli_smoke rfsyn)

  add_test(
    NAME cli_smoke
    COMMAND $<TARGET_FILE:cli_smoke> $<TARGET_FILE:rfsyn>
  )

  set_tests_properties(
    cli_smoke
    PROPERTIES
      LABELS "smoke;cli"
  )

  function(rfsyn_add_validator name group timeout_seconds)
    set(
      validator_source
      "${CMAKE_CURRENT_SOURCE_DIR}/engine/fdtd/${name}.c"
    )

    if(NOT EXISTS "${validator_source}")
      message(FATAL_ERROR "Validator source missing: ${validator_source}")
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
      COMMAND $<TARGET_FILE:${name}>
    )

    set_tests_properties(
      ${name}
      PROPERTIES
        WORKING_DIRECTORY "${validator_workdir}"
        LABELS "fdtd;validation;${group}"
        TIMEOUT "${timeout_seconds}"
    )
  endfunction()

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

  if(RFSYN_BUILD_EXTENDED_VALIDATORS)
    set(
      RFSYN_EXTENDED_VALIDATORS
      validate_cpml_oblique
      validate_materials_high_contrast
      validate_long_term_stability
    )

    # Check all files first so configuration reports every missing source.
    set(RFSYN_MISSING_VALIDATORS)

    foreach(name IN LISTS RFSYN_EXTENDED_VALIDATORS)
      if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/engine/fdtd/${name}.c")
        list(APPEND RFSYN_MISSING_VALIDATORS "engine/fdtd/${name}.c")
      endif()
    endforeach()

    if(RFSYN_MISSING_VALIDATORS)
      string(
        REPLACE ";" "\n  "
        missing_validator_text
        "${RFSYN_MISSING_VALIDATORS}"
      )

      message(
        FATAL_ERROR
        "Extended validators requested, but these files are missing:\n"
        "  ${missing_validator_text}\n"
        "Add the source files or configure with "
        "-DRFSYN_BUILD_EXTENDED_VALIDATORS=OFF."
      )
    endif()

    foreach(name IN LISTS RFSYN_EXTENDED_VALIDATORS)
      rfsyn_add_validator(
        ${name}
        extended
        "${RFSYN_EXTENDED_TIMEOUT}"
      )
    endforeach()
  endif()

elseif(RFSYN_BUILD_EXTENDED_VALIDATORS)
  message(
    FATAL_ERROR
    "RFSYN_BUILD_EXTENDED_VALIDATORS requires BUILD_TESTING=ON"
  )
endif()

message(STATUS "rfsyn: C compiler = ${CMAKE_C_COMPILER}")
message(STATUS "rfsyn: tests = ${BUILD_TESTING}")
message(STATUS "rfsyn: validation refinement = ${RFSYN_VALIDATION_REFINEMENT}")
message(STATUS "rfsyn: extended validators = ${RFSYN_BUILD_EXTENDED_VALIDATORS}")
message(STATUS "rfsyn: sanitizers = ${RFSYN_SANITIZE}")