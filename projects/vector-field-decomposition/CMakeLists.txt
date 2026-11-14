cmake_minimum_required(VERSION 3.10.0)

project(vector-field-decomposition)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Print the build type
# set(CMAKE_BUILD_TYPE "Debug")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()
message(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")

### Configure the compiler
# This is a basic, decent setup that should do something sane on most compilers

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++17 -Wall -Wextra")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-unused-private-field")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3 -fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Setting gcc-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()
  
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  
  # https://gitlab.kitware.com/cmake/cmake/-/issues/22484
  include(CheckCXXCompilerFlag)
  unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
  CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG")
  else()
    SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -mcpu=apple-m1 -DNDEBUG") # Apple M1
  endif()

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)
else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()

add_subdirectory(../../deps/geometry-central deps/geometry-central)
add_subdirectory(../../deps/polyscope deps/polyscope)
add_subdirectory(../../deps/googletest deps/googletest)
set_property(TARGET geometry-central PROPERTY CXX_STANDARD 17)

set(PARAM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../../conformal-parameterization")
set(MODEL_SRC "${CMAKE_CURRENT_SOURCE_DIR}/../principal-curvature-field/src")
set(MESH_CORE
  "${PARAM_ROOT}/QuadMeshing/GlobalFieldsParam/GlobalFieldsParameterization.cpp"
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

# Edit files to build here.
set(SRCS
  src/main.cpp
  src/hodge.cpp
  src/vector_field_lab.cpp
  src/tree-cotree.cpp
  src/harmonic-bases.cpp
  src/hodge-decomposition.cpp
  src/trivial-connections.cpp
  "${MODEL_SRC}/models.cpp"
  ../../utils/src/solvers.cpp
  ${MESH_CORE}
  "${PARAM_ROOT}/VectorFileds/NRosyVectorFields.cpp"
  "${PARAM_ROOT}/VectorFileds/ComplexPolyField/ComplexPolyField.cpp"
)

set(DECOMP_SRCS
  src/tree-cotree.cpp
  src/harmonic-bases.cpp
  src/hodge-decomposition.cpp
  ../../utils/src/solvers.cpp
)

add_executable(main "${SRCS}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main)
target_include_directories(main PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/../Tree-Cotree"
  "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include"
  "${MODEL_SRC}"
  "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args"
  "${PARAM_ROOT}/BaseMesh"
  "${PARAM_ROOT}/Parameterization"
  "${PARAM_ROOT}/QuadMeshing/GlobalFieldsParam"
  "${PARAM_ROOT}/VectorFileds/PrincipalCurvatureFields"
  "${PARAM_ROOT}/VectorFileds"
  "${PARAM_ROOT}/VectorFileds/ComplexPolyField"
  "${PARAM_ROOT}/Solvers"
  "${PARAM_ROOT}/QuadMeshing"
  "/opt/homebrew/include/eigen3"
)
target_link_libraries(main geometry-central polyscope)

add_executable(test-decomp src/test-decomp.cpp ${DECOMP_SRCS})
target_include_directories(test-decomp PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/../Tree-Cotree"
  "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include"
)
target_compile_definitions(test-decomp PRIVATE
  "DECOMP_SOLN_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/include/test-decomp-soln.txt\""
)
target_link_libraries(test-decomp geometry-central gtest)

add_executable(test-connect src/test-connect.cpp src/trivial-connections.cpp ${DECOMP_SRCS})
target_include_directories(test-connect PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/../Tree-Cotree"
  "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include"
)
target_compile_definitions(test-connect PRIVATE
  "CONNECT_SOLN_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/include/test-connect-soln.txt\""
  "CONNECT_SPHERE_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/include/test-connect-sphere.txt\""
)
target_link_libraries(test-connect geometry-central gtest)
