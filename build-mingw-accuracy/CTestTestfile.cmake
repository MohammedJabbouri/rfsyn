# CMake generated Testfile for 
# Source directory: C:/Users/binso/3D Objects/rfsyn
# Build directory: C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(world_unit "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/world_unit.exe")
set_tests_properties(world_unit PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;48;add_test;C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;0;")
add_test(material_catalog_unit "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/material_catalog_unit.exe")
set_tests_properties(material_catalog_unit PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;48;add_test;C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;0;")
add_test(material_config_unit "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/material_config_unit.exe")
set_tests_properties(material_config_unit PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;48;add_test;C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;0;")
add_test(simulate_unit "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/simulate_unit.exe")
set_tests_properties(simulate_unit PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;48;add_test;C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;0;")
add_test(materials_unit "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/materials_unit.exe")
set_tests_properties(materials_unit PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;48;add_test;C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;0;")
add_test(rf_signal_unit "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/rf_signal_unit.exe")
set_tests_properties(rf_signal_unit PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;48;add_test;C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;0;")
add_test(cli_smoke "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/cli_smoke.exe" "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/rfsyn.exe")
set_tests_properties(cli_smoke PROPERTIES  _BACKTRACE_TRIPLES "C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;46;add_test;C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;0;")
add_test(validate_cpml "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/validate_cpml.exe")
set_tests_properties(validate_cpml PROPERTIES  LABELS "fdtd;validation" WORKING_DIRECTORY "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/validation/validate_cpml" _BACKTRACE_TRIPLES "C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;58;add_test;C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;0;")
add_test(validate_cpml_margin "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/validate_cpml_margin.exe")
set_tests_properties(validate_cpml_margin PROPERTIES  LABELS "fdtd;validation" WORKING_DIRECTORY "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/validation/validate_cpml_margin" _BACKTRACE_TRIPLES "C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;58;add_test;C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;0;")
add_test(validate_cpml_planewave "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/validate_cpml_planewave.exe")
set_tests_properties(validate_cpml_planewave PROPERTIES  LABELS "fdtd;validation" WORKING_DIRECTORY "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/validation/validate_cpml_planewave" _BACKTRACE_TRIPLES "C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;58;add_test;C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;0;")
add_test(validate_materials_fresnel "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/validate_materials_fresnel.exe")
set_tests_properties(validate_materials_fresnel PROPERTIES  LABELS "fdtd;validation" WORKING_DIRECTORY "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/validation/validate_materials_fresnel" _BACKTRACE_TRIPLES "C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;58;add_test;C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;0;")
add_test(validate_probe_dft "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/validate_probe_dft.exe")
set_tests_properties(validate_probe_dft PROPERTIES  LABELS "fdtd;validation" WORKING_DIRECTORY "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/validation/validate_probe_dft" _BACKTRACE_TRIPLES "C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;58;add_test;C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;0;")
add_test(validate_materials_energy "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/validate_materials_energy.exe")
set_tests_properties(validate_materials_energy PROPERTIES  LABELS "fdtd;validation" WORKING_DIRECTORY "C:/Users/binso/3D Objects/rfsyn/build-mingw-accuracy/validation/validate_materials_energy" _BACKTRACE_TRIPLES "C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;58;add_test;C:/Users/binso/3D Objects/rfsyn/CMakeLists.txt;0;")
