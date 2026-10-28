cmake_minimum_required(VERSION 3.10.0)

project(principal-curvature-field)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()
message(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++17 -Wall -Wextra")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-unused-private-field")
  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG "-g3 -fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    SET(CMAKE_CXX_FLAGS_DEBUG "-g3")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  include(CheckCXXCompilerFlag)
  unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
  CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
  else()
    SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -mcpu=apple-m1 -DNDEBUG")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  set(CMAKE_CXX_FLAGS "/W4 /MP /std:c++17")
  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)
else()
  message(FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]")
endif()

add_subdirectory(../../deps/geometry-central deps/geometry-central)
add_subdirectory(../../deps/polyscope deps/polyscope)
# This geometry-central pins C++11, but the installed Eigen needs C++14.
set_property(TARGET geometry-central PROPERTY CXX_STANDARD 17)

set(PARAM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../../conformal-parameterization")
set(MESH_SOURCES
  "${PARAM_ROOT}/Parameterization/GlobalFieldsParam/QuadCover.cpp"
  "${PARAM_ROOT}/Parameterization/GlobalFieldsParam/GlobalFieldsParameterization.cpp"
  "${PARAM_ROOT}/Parameterization/Parameterization.cpp"
  "${PARAM_ROOT}/VectorFileds/PrincipalCurvatureFields/PrincipalCurvatureField.cpp"
  "${PARAM_ROOT}/BaseMesh/Mesh.cpp"
  "${PARAM_ROOT}/BaseMesh/MeshIO.cpp"
  "${PARAM_ROOT}/BaseMesh/Vertex.cpp"
  "${PARAM_ROOT}/BaseMesh/Edge.cpp"
  "${PARAM_ROOT}/BaseMesh/Face.cpp"
  "${PARAM_ROOT}/BaseMesh/HalfEdge.cpp"
  "${PARAM_ROOT}/BaseMesh/QcError.cpp"
)

function(link_viewer target)
  target_include_directories(${target} PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args"
    "${PARAM_ROOT}/BaseMesh"
    "${PARAM_ROOT}/Parameterization"
    "${PARAM_ROOT}/Parameterization/GlobalFieldsParam"
    "${PARAM_ROOT}/VectorFileds/PrincipalCurvatureFields"
    "${PARAM_ROOT}/VectorFileds"
    "${PARAM_ROOT}/VectorFileds/ComplexPolyField"
    "${PARAM_ROOT}/Solvers"
    "${PARAM_ROOT}/QuadMeshing"
    "/opt/homebrew/include/eigen3"
  )
  target_link_libraries(${target} geometry-central polyscope)
endfunction()

add_executable(main
  src/main.cpp
  src/models.cpp
  ${MESH_SOURCES}
  "${PARAM_ROOT}/VectorFileds/NRosyVectorFields.cpp"
  "${PARAM_ROOT}/VectorFileds/ComplexPolyField/ComplexPolyField.cpp"
  "${PARAM_ROOT}/Parameterization/GlobalFieldsParam/MIQQuad.cpp"
  "${PARAM_ROOT}/Solvers/MixedIntegerProgram.cpp"
)
link_viewer(main)

add_executable(quad_test
  src/quad_test.cpp
  src/models.cpp
  ${MESH_SOURCES}
  "${PARAM_ROOT}/Parameterization/GlobalFieldsParam/MIQQuad.cpp"
  "${PARAM_ROOT}/Solvers/MixedIntegerProgram.cpp"
  "${PARAM_ROOT}/QuadMeshing/NaiveQuadMerge.cpp"
)
link_viewer(quad_test)

enable_language(OBJCXX)
set(CMAKE_OBJCXX_STANDARD 17)

add_executable(simple_param
  ../simple-param/src/main.mm
  "${PARAM_ROOT}/Parameterization/Parameterization.cpp"
  "${PARAM_ROOT}/Parameterization/SimpleParam/LSCM/Lscm.cpp"
  "${PARAM_ROOT}/Parameterization/SimpleParam/LSCM/Scp.cpp"
  "${PARAM_ROOT}/Parameterization/SimpleParam/ARAP/ARAP.cpp"
  "${PARAM_ROOT}/Parameterization/SimpleParam/ARAP/Tutte.cpp"
  "${PARAM_ROOT}/Parameterization/SimpleParam/ABF/AbfPlusPlus.cpp"
  "${PARAM_ROOT}/Parameterization/SimpleParam/ABF/LinAbf.cpp"
  "${PARAM_ROOT}/Solvers/AugmentedLagrangian.cpp"
  "${PARAM_ROOT}/Bounded/bounded_distortion_mapping/BoundedDistortionMapping.cpp"
  "${PARAM_ROOT}/BaseMesh/Mesh.cpp"
  "${PARAM_ROOT}/BaseMesh/MeshIO.cpp"
  "${PARAM_ROOT}/BaseMesh/Vertex.cpp"
  "${PARAM_ROOT}/BaseMesh/Edge.cpp"
  "${PARAM_ROOT}/BaseMesh/Face.cpp"
  "${PARAM_ROOT}/BaseMesh/HalfEdge.cpp"
  "${PARAM_ROOT}/BaseMesh/QcError.cpp"
  "${PARAM_ROOT}/CutSeamMesh/CutSeamMesh.cpp"
)
target_include_directories(simple_param PRIVATE
  "${PARAM_ROOT}/BaseMesh"
  "${PARAM_ROOT}/Parameterization"
  "${PARAM_ROOT}/Solvers"
  "${PARAM_ROOT}/Parameterization/SimpleParam/LSCM"
  "${PARAM_ROOT}/Parameterization/SimpleParam/ARAP"
  "${PARAM_ROOT}/Parameterization/SimpleParam/ABF"
  "${PARAM_ROOT}/Bounded/bounded_distortion_mapping"
  "${PARAM_ROOT}/CutSeamMesh"
  "/opt/homebrew/include/eigen3"
)
target_compile_options(simple_param PRIVATE "-fobjc-arc" "-std=c++17")
target_link_libraries(simple_param glfw "-framework OpenGL" "-framework Cocoa" "-framework CoreText")

set(OSQP_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/osqp")
set(UNITTESTS OFF CACHE BOOL "OSQP unit tests" FORCE)
set(ENABLE_MKL_PARDISO OFF CACHE BOOL "OSQP MKL Pardiso" FORCE)
add_subdirectory(${OSQP_ROOT} deps/osqp)

add_executable(cone_unwrap
  ../cone-unwrap/src/main.mm
  src/models.cpp
  "${PARAM_ROOT}/Parameterization/Parameterization.cpp"
  "${PARAM_ROOT}/Parameterization/ConeParam/ConeParameterization.cpp"
  "${PARAM_ROOT}/Parameterization/ConeParam/Conformal/CirclePatterns.cpp"
  "${PARAM_ROOT}/Parameterization/ConeParam/Conformal/Cetm.cpp"
  "${PARAM_ROOT}/Parameterization/ConeParam/RicciFlow/RicciFlow.cpp"
  "${PARAM_ROOT}/Parameterization/ConeParam/IncrementalFlattenning/IncrementalFlattening.cpp"
  "${PARAM_ROOT}/Parameterization/SimpleParam/ARAP/Tutte.cpp"
  "${PARAM_ROOT}/Solvers/Solver.cpp"
  "${PARAM_ROOT}/BaseMesh/Mesh.cpp"
  "${PARAM_ROOT}/BaseMesh/MeshIO.cpp"
  "${PARAM_ROOT}/BaseMesh/Vertex.cpp"
  "${PARAM_ROOT}/BaseMesh/Edge.cpp"
  "${PARAM_ROOT}/BaseMesh/Face.cpp"
  "${PARAM_ROOT}/BaseMesh/HalfEdge.cpp"
  "${PARAM_ROOT}/BaseMesh/QcError.cpp"
  "${PARAM_ROOT}/BaseMesh/GaussianCurvature.cpp"
)
target_include_directories(cone_unwrap PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${PARAM_ROOT}/BaseMesh"
  "${PARAM_ROOT}/Parameterization"
  "${PARAM_ROOT}/Parameterization/ConeParam"
  "${PARAM_ROOT}/Parameterization/ConeParam/Conformal"
  "${PARAM_ROOT}/Parameterization/ConeParam/RicciFlow"
  "${PARAM_ROOT}/Parameterization/ConeParam/IncrementalFlattenning"
  "${PARAM_ROOT}/Parameterization/SimpleParam/ARAP"
  "${PARAM_ROOT}/CutSeamMesh"
  "${PARAM_ROOT}/Solvers"
  "${PARAM_ROOT}/Solvers/Mosek"
  "/opt/homebrew/include/eigen3"
)
target_compile_definitions(cone_unwrap PRIVATE CONE_USE_OSQP)
set_source_files_properties(../cone-unwrap/src/main.mm PROPERTIES COMPILE_FLAGS "-fobjc-arc")
target_link_libraries(cone_unwrap osqpstatic glfw "-framework OpenGL" "-framework Cocoa" "-framework CoreText")

add_executable(hof_unwrap
  ../hof-unwrap/src/main.mm
  src/models.cpp
  "${PARAM_ROOT}/Parameterization/Parameterization.cpp"
  "${PARAM_ROOT}/Parameterization/CutSeamParam/CutSeamParameterization.cpp"
  "${PARAM_ROOT}/Parameterization/CutSeamParam/HoloOneForm/HolomorphicOneForm.cpp"
  "${PARAM_ROOT}/CutSeamMesh/CutSeamMesh.cpp"
  "${PARAM_ROOT}/BaseMesh/Mesh.cpp"
  "${PARAM_ROOT}/BaseMesh/MeshIO.cpp"
  "${PARAM_ROOT}/BaseMesh/Vertex.cpp"
  "${PARAM_ROOT}/BaseMesh/Edge.cpp"
  "${PARAM_ROOT}/BaseMesh/Face.cpp"
  "${PARAM_ROOT}/BaseMesh/HalfEdge.cpp"
  "${PARAM_ROOT}/BaseMesh/QcError.cpp"
)
target_include_directories(hof_unwrap PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${PARAM_ROOT}/BaseMesh"
  "${PARAM_ROOT}/Parameterization"
  "${PARAM_ROOT}/Parameterization/CutSeamParam"
  "${PARAM_ROOT}/Parameterization/CutSeamParam/HoloOneForm"
  "${PARAM_ROOT}/CutSeamMesh"
  "/opt/homebrew/include/eigen3"
)
set_source_files_properties(../hof-unwrap/src/main.mm PROPERTIES COMPILE_FLAGS "-fobjc-arc")
target_link_libraries(hof_unwrap glfw "-framework OpenGL" "-framework Cocoa" "-framework CoreText")
