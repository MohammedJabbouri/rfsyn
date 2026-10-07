cmake_minimum_required(VERSION 3.16)
project(rfsyn C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
  add_compile_options(/W3)
else()
  add_compile_options(-Wall -Wextra)
endif()

option(RFSYN_SANITIZE "Build with -fsanitize=address,undefined (gcc/clang only)" OFF)

find_package(OpenMP QUIET)
if(OpenMP_C_FOUND)
  message(STATUS "rfsyn: OpenMP found -- parallel loops enabled")
else()
  message(STATUS "rfsyn: OpenMP not found -- building single-threaded")
endif()

add_library(rfsyn_fdtd STATIC
  engine/fdtd/grid.c
  engine/fdtd/cpml.c
  engine/fdtd/materials.c
  engine/fdtd/source.c
  engine/fdtd/probe.c
  engine/fdtd/simulate.c
)

target_link_libraries(rfsyn_fdtd PUBLIC m)
if(OpenMP_C_FOUND)
  target_link_libraries(rfsyn_fdtd PUBLIC OpenMP::OpenMP_C)
endif()

add_library(rfsyn_core STATIC
  engine/core/chain.c
  engine/core/chains.c
  engine/core/rf_signal.c
  engine/core/rng.c
)

add_library(rfsyn_world STATIC
  engine/world/world.c
  engine/world/material_catalog.c

)

add_executable(world_unit tests/world_unit.c)
target_link_libraries(world_unit PRIVATE rfsyn_world m)

target_link_libraries(rfsyn_world PUBLIC m)
target_link_libraries(rfsyn_core PUBLIC m)

add_executable(material_catalog_unit tests/material_catalog_unit.c)
target_link_libraries(material_catalog_unit PRIVATE rfsyn_world m)

add_library(rfsyn_noise STATIC
  engine/noise/awgn.c
)
target_link_libraries(rfsyn_noise PUBLIC rfsyn_core m)

add_executable(rfsyn
  cli/main.c
  cli/cJSON.c
  cli/config.c
  cli/presets.c
  cli/settings.c
  cli/fdtd_stage.c
)
target_link_libraries(rfsyn PRIVATE rfsyn_core rfsyn_noise rfsyn_fdtd m)

add_executable(cli_smoke tests/cli_smoke.c)

add_executable(rf_signal_unit tests/rf_signal_unit.c)
target_link_libraries(rf_signal_unit PRIVATE rfsyn_core m)

add_executable(simulate_unit tests/simulate_unit.c)
target_link_libraries(simulate_unit PRIVATE rfsyn_fdtd m)

add_executable(materials_unit tests/materials_unit.c)
target_link_libraries(materials_unit PRIVATE rfsyn_fdtd m)

# FDTD validators 
set(RFSYN_FDTD_VALIDATORS
  validate_cpml
  validate_cpml_margin
  validate_materials_fresnel
  validate_probe_dft
  validate_materials_energy
)

foreach(validator IN LISTS RFSYN_FDTD_VALIDATORS)
  add_executable(${validator} engine/fdtd/${validator}.c)
  target_link_libraries(${validator} PRIVATE rfsyn_fdtd m)
endforeach()

enable_testing()

add_test(NAME simulate_unit COMMAND simulate_unit)
add_test(NAME materials_unit COMMAND materials_unit)
add_test(NAME world_unit COMMAND world_unit)
add_test(NAME material_catalog_unit COMMAND material_catalog_unit)

foreach(test_name IN ITEMS validate_cpml_margin validate_materials_fresnel validate_probe_dft validate_materials_energy)
  add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()

add_test(NAME validate_cpml_smoke COMMAND validate_cpml)

add_test(NAME cli_smoke COMMAND cli_smoke $<TARGET_FILE:rfsyn>)

add_test(NAME rf_signal_unit COMMAND rf_signal_unit)

if(RFSYN_SANITIZE)
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    set(RFSYN_ALL_TARGETS rfsyn_fdtd rfsyn_core rfsyn_noise rfsyn_world rfsyn cli_smoke rf_signal_unit simulate_unit materials_unit world_unit ${RFSYN_FDTD_VALIDATORS})
    foreach(tgt IN LISTS RFSYN_ALL_TARGETS)
      target_compile_options(${tgt} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer -g)
      target_link_options(${tgt} PRIVATE -fsanitize=address,undefined)
    endforeach()
  else()
    message(WARNING "rfsyn: RFSYN_SANITIZE requested but compiler is not GNU/Clang -- ignoring")
  endif()
endif()