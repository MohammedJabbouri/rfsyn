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
set_property(CACHE RFSYN_VALIDATION_REFINEMENT PROPERTY STRINGS 1 2 4)

set(RFSYN_STABILITY_STEPS "100000" CACHE STRING "Steps per stability case")
set(RFSYN_VALIDATION_TIMEOUT "900" CACHE STRING "Baseline validator timeout")
set(RFSYN_EXTENDED_TIMEOUT "7200" CACHE STRING "Extended validator timeout")

if(NOT RFSYN_VALIDATION_REFINEMENT MATCHES "^(1|2|4)$")
  message(FATAL_ERROR "RFSYN_VALIDATION_REFINEMENT must be 1, 2, or 4")
endif()

if(RFSYN_STABILITY_STEPS LESS 10000 OR RFSYN_STABILITY_STEPS GREATER 10000000)
  message(FATAL_ERROR "RFSYN_STABILITY_STEPS must be between 10000 and 10000000")
endif()

foreach(timeout IN ITEMS RFSYN_VALIDATION_TIMEOUT RFSYN_EXTENDED_TIMEOUT)
  if(NOT ${timeout} MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "${timeout} must be a positive integer")
  endif()
endforeach()

add_library(rfsyn_build_options INTERFACE)

if(RFSYN_SANITIZE)
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang" AND NOT MSVC)
    target_compile_options(rfsyn_build_options INTERFACE
      -fsanitize=address,undefined -fno-omit-frame-pointer -g
    )
    target_link_options(rfsyn_build_options INTERFACE
      -fsanitize=address,undefined
    )
  else()
    message(FATAL_ERROR "Sanitizers require GCC or Clang")
  endif()
endif()

set(RFSYN_MATH_LIBRARIES)
if(NOT MSVC)
  list(APPEND RFSYN_MATH_LIBRARIES m)
endif()

function(rfsyn_apply_build_options target)
  target_link_libraries(${target} PRIVATE rfsyn_build_options)
endfunction()

add_library(rfsyn_fdtd STATIC
  engine/fdtd/grid.c
  engine/fdtd/cpml.c
  engine/fdtd/materials.c
  engine/fdtd/source.c
  engine/fdtd/probe.c
  engine/fdtd/simulate.c
  engine/fdtd/subpixel_te.c
)

add_library(rfsyn_core STATIC
  engine/core/chain.c
  engine/core/chains.c
  engine/core/rf_signal.c
  engine/core/rng.c
)

add_library(rfsyn_cjson STATIC cli/cJSON.c)
target_include_directories(rfsyn_cjson PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/cli")

add_library(rfsyn_world STATIC
  engine/world/world.c
  engine/world/material_catalog.c
  cli/material_config.c
)
target_include_directories(rfsyn_world PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

add_library(rfsyn_noise STATIC engine/noise/awgn.c)

target_link_libraries(rfsyn_fdtd PUBLIC ${RFSYN_MATH_LIBRARIES})
target_link_libraries(rfsyn_core PUBLIC ${RFSYN_MATH_LIBRARIES})
target_link_libraries(rfsyn_world PUBLIC rfsyn_cjson ${RFSYN_MATH_LIBRARIES})
target_link_libraries(rfsyn_noise PUBLIC rfsyn_core ${RFSYN_MATH_LIBRARIES})

foreach(lib IN ITEMS rfsyn_fdtd rfsyn_core rfsyn_cjson rfsyn_world rfsyn_noise)
  rfsyn_apply_build_options(${lib})
endforeach()

if(RFSYN_ENABLE_OPENMP)
  find_package(OpenMP QUIET COMPONENTS C)
  if(OpenMP_C_FOUND)
    target_link_libraries(rfsyn_fdtd PUBLIC OpenMP::OpenMP_C)
    target_compile_definitions(rfsyn_fdtd PUBLIC RFSYN_HAS_OPENMP=1)
  endif()
endif()

add_executable(rfsyn
  cli/main.c
  cli/settings.c
  cli/config.c
  cli/presets.c
  cli/fdtd_stage.c
)
target_link_libraries(rfsyn PRIVATE
  rfsyn_fdtd rfsyn_core rfsyn_noise rfsyn_world ${RFSYN_MATH_LIBRARIES}
)
rfsyn_apply_build_options(rfsyn)

if(RFSYN_BUILD_PERIODIC_ENERGY_AUDIT)
  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/engine/fdtd/audit_periodic_energy.c")
    message(FATAL_ERROR
      "Missing engine/fdtd/audit_periodic_energy.c. "
      "Set -DRFSYN_BUILD_PERIODIC_ENERGY_AUDIT=OFF."
    )
  endif()
  add_executable(audit_periodic_energy engine/fdtd/audit_periodic_energy.c)
  target_link_libraries(audit_periodic_energy PRIVATE
    rfsyn_fdtd ${RFSYN_MATH_LIBRARIES}
  )
  rfsyn_apply_build_options(audit_periodic_energy)
endif()

include(CTest)

if(BUILD_TESTING)

  function(rfsyn_add_unit_test name)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/${name}.c")
      message(FATAL_ERROR "Missing unit-test source: tests/${name}.c")
    endif()
    add_executable(${name} "tests/${name}.c")
    target_link_libraries(${name} PRIVATE ${ARGN} ${RFSYN_MATH_LIBRARIES})
    rfsyn_apply_build_options(${name})
    add_test(NAME ${name} COMMAND $<TARGET_FILE:${name}>)
    set_tests_properties(${name} PROPERTIES LABELS "unit")
  endfunction()

  function(rfsyn_add_validator name group timeout_seconds)
    set(src "${CMAKE_CURRENT_SOURCE_DIR}/engine/fdtd/${name}.c")
    if(NOT EXISTS "${src}")
      message(FATAL_ERROR "Missing validator source: ${src}")
    endif()
    add_executable(${name} "${src}")
    target_link_libraries(${name} PRIVATE rfsyn_fdtd ${RFSYN_MATH_LIBRARIES})
    target_compile_definitions(${name} PRIVATE
      RFSYN_VALIDATION_REFINEMENT=${RFSYN_VALIDATION_REFINEMENT}
    )
    rfsyn_apply_build_options(${name})
    set(workdir "${CMAKE_CURRENT_BINARY_DIR}/validation/${name}")
    file(MAKE_DIRECTORY "${workdir}")
    add_test(NAME ${name} COMMAND $<TARGET_FILE:${name}> ${ARGN})
    if(group STREQUAL "measurement")
      set(labels "fdtd;measurement;tilted")
    else()
      set(labels "fdtd;validation;${group}")
    endif()
    set_tests_properties(${name} PROPERTIES
      WORKING_DIRECTORY "${workdir}"
      LABELS "${labels}"
      TIMEOUT "${timeout_seconds}"
    )
  endfunction()

  rfsyn_add_unit_test(world_unit rfsyn_world)
  rfsyn_add_unit_test(material_catalog_unit rfsyn_world)
  rfsyn_add_unit_test(material_config_unit rfsyn_world)
  rfsyn_add_unit_test(simulate_unit rfsyn_fdtd)
  rfsyn_add_unit_test(materials_unit rfsyn_fdtd)
  rfsyn_add_unit_test(rf_signal_unit rfsyn_core)
  rfsyn_add_unit_test(subpixel_te_unit rfsyn_fdtd)
  rfsyn_add_unit_test(periodic_materials_unit rfsyn_fdtd)

  set_tests_properties(periodic_materials_unit PROPERTIES LABELS "unit;fdtd;periodic")

  add_executable(cli_smoke tests/cli_smoke.c)
  target_link_libraries(cli_smoke PRIVATE
    rfsyn_core rfsyn_noise rfsyn_fdtd rfsyn_world ${RFSYN_MATH_LIBRARIES}
  )
  rfsyn_apply_build_options(cli_smoke)
  add_dependencies(cli_smoke rfsyn)
  add_test(NAME cli_smoke COMMAND $<TARGET_FILE:cli_smoke> $<TARGET_FILE:rfsyn>)
  set_tests_properties(cli_smoke PROPERTIES LABELS "smoke;cli")

  foreach(v IN ITEMS
      validate_cpml
      validate_cpml_margin
      validate_cpml_planewave
      validate_materials_fresnel
      validate_probe_dft
      validate_materials_energy)
    rfsyn_add_validator(${v} baseline "${RFSYN_VALIDATION_TIMEOUT}")
  endforeach()

  if(RFSYN_BUILD_EXTENDED_VALIDATORS)
    rfsyn_add_validator(validate_cpml_oblique extended "${RFSYN_EXTENDED_TIMEOUT}")
    rfsyn_add_validator(validate_materials_high_contrast extended "${RFSYN_EXTENDED_TIMEOUT}" --full)
    rfsyn_add_validator(validate_long_term_stability extended "${RFSYN_EXTENDED_TIMEOUT}" --steps "${RFSYN_STABILITY_STEPS}")
    rfsyn_add_validator(validate_materials_tilted measurement "${RFSYN_EXTENDED_TIMEOUT}")
  endif()

endif()

message(STATUS "rfsyn: C compiler = ${CMAKE_C_COMPILER}")
message(STATUS "rfsyn: tests = ${BUILD_TESTING}")
message(STATUS "rfsyn: refinement = ${RFSYN_VALIDATION_REFINEMENT}")
message(STATUS "rfsyn: OpenMP option = ${RFSYN_ENABLE_OPENMP}")
message(STATUS "rfsyn: extended validators = ${RFSYN_BUILD_EXTENDED_VALIDATORS}")
message(STATUS "rfsyn: periodic-energy audit = ${RFSYN_BUILD_PERIODIC_ENERGY_AUDIT}")
message(STATUS "rfsyn: stability steps = ${RFSYN_STABILITY_STEPS}")
message(STATUS "rfsyn: sanitizers = ${RFSYN_SANITIZE}")