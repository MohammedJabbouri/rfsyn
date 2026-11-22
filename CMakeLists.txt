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
option(RFSYN_BUILD_EXTENDED_VALIDATORS "Build extended validators" ON)
option(RFSYN_BUILD_PERIODIC_ENERGY_AUDIT "Build periodic-energy audit" ON)

set(RFSYN_VALIDATION_REFINEMENT "1" CACHE STRING "Validator refinement: 1, 2, or 4")
set(RFSYN_STABILITY_STEPS "1000" CACHE STRING "Steps per stability case")
set(RFSYN_VALIDATION_TIMEOUT "900" CACHE STRING "Baseline validator timeout")
set(RFSYN_EXTENDED_TIMEOUT "7200" CACHE STRING "Extended validator timeout")

set_property(CACHE RFSYN_VALIDATION_REFINEMENT PROPERTY STRINGS 1 2 4)

if(NOT "${RFSYN_VALIDATION_REFINEMENT}" MATCHES "^(1|2|4)$")
  message(FATAL_ERROR "RFSYN_VALIDATION_REFINEMENT must be 1, 2, or 4")
endif()

foreach(variable IN ITEMS
        RFSYN_STABILITY_STEPS
        RFSYN_VALIDATION_TIMEOUT
        RFSYN_EXTENDED_TIMEOUT)
  if(NOT "${${variable}}" MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "${variable} must be a positive integer")
  endif()
endforeach()

if(RFSYN_STABILITY_STEPS LESS 10000 OR
        RFSYN_STABILITY_STEPS GREATER 10000000)
  message(FATAL_ERROR "RFSYN_STABILITY_STEPS must be between 10000 and 10000000")
endif()

include(CTest)

add_library(rfsyn_build_options INTERFACE
        tests/materials_tensor_unit.c)

if(RFSYN_SANITIZE)
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang" AND NOT MSVC)
    target_compile_options(rfsyn_build_options INTERFACE
            -fsanitize=address,undefined
            -fno-omit-frame-pointer
            -g
    )
    target_link_options(rfsyn_build_options INTERFACE
            -fsanitize=address,undefined
    )
  else()
    message(FATAL_ERROR "Sanitizers require a supported GCC or Clang toolchain")
  endif()
endif()

set(RFSYN_MATH_LIBRARIES)

if(NOT MSVC)
  list(APPEND RFSYN_MATH_LIBRARIES m)
endif()

function(rfsyn_require_files)
  foreach(source IN LISTS ARGN)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
      message(FATAL_ERROR "Missing required source file: ${source}")
    endif()
  endforeach()
endfunction()

function(rfsyn_add_library name)
  rfsyn_require_files(${ARGN})
  add_library(${name} STATIC ${ARGN})
  target_link_libraries(${name} PRIVATE rfsyn_build_options)
endfunction()

function(rfsyn_add_executable name)
  rfsyn_require_files(${ARGN})
  add_executable(${name} ${ARGN})
  target_link_libraries(${name} PRIVATE rfsyn_build_options)
endfunction()

rfsyn_add_library(rfsyn_fdtd
        engine/fdtd/grid.c
        engine/fdtd/cpml.c
        engine/fdtd/materials.c
        engine/fdtd/source.c
        engine/fdtd/probe.c
        engine/fdtd/simulate.c
        engine/fdtd/subpixel_te.c
)

rfsyn_add_library(rfsyn_core
        engine/core/chain.c
        engine/core/chains.c
        engine/core/rf_signal.c
        engine/core/rng.c
)

rfsyn_add_library(rfsyn_cjson
        cli/cJSON.c
)

rfsyn_add_library(rfsyn_world
        engine/world/world.c
        engine/world/material_catalog.c
        cli/material_config.c
)

rfsyn_add_library(rfsyn_noise
        engine/noise/awgn.c
)

target_include_directories(rfsyn_cjson PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/cli"
)

target_include_directories(rfsyn_world PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

if(RFSYN_MATH_LIBRARIES)
  foreach(target IN ITEMS rfsyn_fdtd rfsyn_core rfsyn_world rfsyn_noise)
    target_link_libraries(${target} PUBLIC ${RFSYN_MATH_LIBRARIES})
  endforeach()
endif()

target_link_libraries(rfsyn_world PUBLIC rfsyn_cjson)
target_link_libraries(rfsyn_noise PUBLIC rfsyn_core)

set(RFSYN_OPENMP_ACTIVE OFF)

if(RFSYN_ENABLE_OPENMP)
  find_package(OpenMP QUIET COMPONENTS C)

  if(OpenMP_C_FOUND)
    target_link_libraries(rfsyn_fdtd PUBLIC OpenMP::OpenMP_C)
    target_compile_definitions(rfsyn_fdtd PUBLIC RFSYN_HAS_OPENMP=1)
    set(RFSYN_OPENMP_ACTIVE ON)
  endif()
endif()

rfsyn_add_executable(rfsyn
        cli/main.c
        cli/settings.c
        cli/config.c
        cli/presets.c
        cli/fdtd_stage.c
)

target_link_libraries(rfsyn PRIVATE
        rfsyn_fdtd
        rfsyn_core
        rfsyn_noise
        rfsyn_world
        ${RFSYN_MATH_LIBRARIES}
)

if(RFSYN_BUILD_PERIODIC_ENERGY_AUDIT)
  rfsyn_require_files(
          engine/fdtd/audit_periodic_energy.c
          engine/fdtd/extended_test.h
  )

  rfsyn_add_executable(audit_periodic_energy
          engine/fdtd/audit_periodic_energy.c
  )

  target_link_libraries(audit_periodic_energy PRIVATE
          rfsyn_fdtd
          ${RFSYN_MATH_LIBRARIES}
  )
endif()

if(BUILD_TESTING)
  function(rfsyn_add_unit_test name)
    rfsyn_add_executable(${name} "tests/${name}.c")

    target_link_libraries(${name} PRIVATE
            ${ARGN}
            ${RFSYN_MATH_LIBRARIES}
    )

    add_test(NAME ${name} COMMAND $<TARGET_FILE:${name}>)
    set_tests_properties(${name} PROPERTIES LABELS "unit")
  endfunction()

  function(rfsyn_add_validator name group timeout_seconds)
    rfsyn_add_executable(${name} "engine/fdtd/${name}.c")

    target_link_libraries(${name} PRIVATE
            rfsyn_fdtd
            ${RFSYN_MATH_LIBRARIES}
    )

    target_compile_definitions(${name} PRIVATE
            RFSYN_VALIDATION_REFINEMENT=${RFSYN_VALIDATION_REFINEMENT}
    )

    set(workdir "${CMAKE_CURRENT_BINARY_DIR}/validation/${name}")
    file(MAKE_DIRECTORY "${workdir}")

    if("${group}" STREQUAL "measurement")
      set(labels "fdtd;measurement;tilted")
    else()
      set(labels "fdtd;validation;${group}")
    endif()

    add_test(NAME ${name} COMMAND $<TARGET_FILE:${name}> ${ARGN})

    set_tests_properties(${name} PROPERTIES
            WORKING_DIRECTORY "${workdir}"
            LABELS "${labels}"
            TIMEOUT "${timeout_seconds}"
    )
  endfunction()

  foreach(name IN ITEMS
          world_unit
          material_catalog_unit
          material_config_unit)
    rfsyn_add_unit_test(${name} rfsyn_world)
  endforeach()

  foreach(name IN ITEMS
          simulate_unit
          materials_unit
          materials_tensor_unit
          subpixel_te_unit
          periodic_materials_unit
          periodic_wrapping_unit)
    rfsyn_add_unit_test(${name} rfsyn_fdtd)
  endforeach()

  rfsyn_add_unit_test(rf_signal_unit rfsyn_core)

  set_tests_properties(
          periodic_materials_unit
          periodic_wrapping_unit
          PROPERTIES LABELS "unit;fdtd;periodic"
  )

  rfsyn_add_executable(cli_smoke tests/cli_smoke.c)

  target_link_libraries(cli_smoke PRIVATE
          rfsyn_core
          rfsyn_noise
          rfsyn_fdtd
          rfsyn_world
          ${RFSYN_MATH_LIBRARIES}
  )

  add_dependencies(cli_smoke rfsyn)

  add_test(NAME cli_smoke COMMAND
          $<TARGET_FILE:cli_smoke>
          $<TARGET_FILE:rfsyn>
  )

  set_tests_properties(cli_smoke PROPERTIES LABELS "smoke;cli")

  foreach(name IN ITEMS
          validate_cpml
          validate_cpml_margin
          validate_cpml_planewave
          validate_materials_fresnel
          validate_probe_dft
          validate_materials_energy)
    rfsyn_add_validator(${name} baseline "${RFSYN_VALIDATION_TIMEOUT}")
  endforeach()

  if(RFSYN_BUILD_EXTENDED_VALIDATORS)
    rfsyn_require_files(engine/fdtd/extended_test.h)

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
            --steps "${RFSYN_STABILITY_STEPS}"
    )

    rfsyn_add_validator(
            validate_materials_tilted
            measurement
            "${RFSYN_EXTENDED_TIMEOUT}"
    )
  endif()
endif()

message(STATUS "rfsyn: C compiler = ${CMAKE_C_COMPILER}")
message(STATUS "rfsyn: tests = ${BUILD_TESTING}")
message(STATUS "rfsyn: refinement = ${RFSYN_VALIDATION_REFINEMENT}")
message(STATUS "rfsyn: OpenMP requested = ${RFSYN_ENABLE_OPENMP}")
message(STATUS "rfsyn: OpenMP active = ${RFSYN_OPENMP_ACTIVE}")
message(STATUS "rfsyn: extended validators option = ${RFSYN_BUILD_EXTENDED_VALIDATORS}")
message(STATUS "rfsyn: periodic-energy audit = ${RFSYN_BUILD_PERIODIC_ENERGY_AUDIT}")
message(STATUS "rfsyn: stability steps = ${RFSYN_STABILITY_STEPS}")
message(STATUS "rfsyn: sanitizers = ${RFSYN_SANITIZE}")