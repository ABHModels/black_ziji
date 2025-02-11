cmake_minimum_required(VERSION 3.10)
project(ZijirayProject VERSION 1.0 LANGUAGES CXX)

# Force shared libraries to have a .so extension on macOS.
if(APPLE)
  set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
endif()

# Allow the user to force using Homebrew libraries on macOS.
if(APPLE)
  option(USE_HOMEBREW_LIBS "Use Homebrew-installed libraries on macOS" ON)
else()
  option(USE_HOMEBREW_LIBS "Use Homebrew-installed libraries on macOS" OFF)
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Platform-specific configuration ---
if(APPLE)
  message(STATUS "Platform: macOS")

  if(USE_HOMEBREW_LIBS)
    # Get Homebrew prefix for libomp
    execute_process(
      COMMAND brew --prefix libomp
      OUTPUT_VARIABLE HB_LIBOMP_PREFIX
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    message(STATUS "Homebrew libomp prefix: ${HB_LIBOMP_PREFIX}")

    # Get Homebrew prefix for GSL
    execute_process(
      COMMAND brew --prefix gsl
      OUTPUT_VARIABLE HB_GSL_PREFIX
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    message(STATUS "Homebrew GSL prefix: ${HB_GSL_PREFIX}")

    # Optionally, force using a Homebrew-installed g++.
    if(NOT DEFINED CMAKE_CXX_COMPILER_USER)
      # Adjust the compiler name if needed (e.g. g++-12, g++-11, etc.)
      set(CMAKE_CXX_COMPILER "/opt/homebrew/bin/g++-12")
      message(STATUS "Using Homebrew's g++: ${CMAKE_CXX_COMPILER}")
    else()
      message(STATUS "Using user-provided C++ compiler: ${CMAKE_CXX_COMPILER}")
    endif()
    
    # Add include and link directories from Homebrew
    if(HB_LIBOMP_PREFIX)
      include_directories("${HB_LIBOMP_PREFIX}/include")
      link_directories("${HB_LIBOMP_PREFIX}/lib")
    else()
      message(WARNING "Could not determine Homebrew prefix for libomp. Check your Homebrew installation.")
    endif()

    if(HB_GSL_PREFIX)
      include_directories("${HB_GSL_PREFIX}/include")
      link_directories("${HB_GSL_PREFIX}/lib")
    else()
      message(WARNING "Could not determine Homebrew prefix for gsl. Check your Homebrew installation.")
    endif()

  else()
    message(STATUS "Using system libraries on macOS (not Homebrew)")
  endif()
else()
  message(STATUS "Platform: Linux/Other")
endif()

# --- Find Packages (OpenMP, GSL) ---
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARY_DIRS})

# --- Set OpenMP flags ---
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -fopenmp")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
endif()

# --- Include Project Directories ---
include_directories(
  ${CMAKE_SOURCE_DIR}/zijiray/include
  ${CMAKE_SOURCE_DIR}/external/xtensor/include
  ${CMAKE_SOURCE_DIR}/external/xtl/include
)

# --- Define Source Files ---
file(GLOB SRC_RAY
  "${CMAKE_SOURCE_DIR}/zijiray/src/utils/*.cc"
  "${CMAKE_SOURCE_DIR}/zijiray/src/GR/*.cc"
  "${CMAKE_SOURCE_DIR}/zijiray/src/raytrace/*.cc"
)

set(TRANSIT_SRC ${SRC_RAY}
  ${CMAKE_SOURCE_DIR}/zijiray/src/transit/main.cc
)
set(LPG_SRC ${SRC_RAY}
  ${CMAKE_SOURCE_DIR}/zijiray/src/lpg/main.cc
)

# --- Set Output Directory for Shared Libraries ---
set(LIB_OUTPUT_DIR ${CMAKE_BINARY_DIR}/bin/ray)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_OUTPUT_DIR})

# --- Build Targets ---
add_library(transit_ray SHARED ${TRANSIT_SRC})
add_library(lpg SHARED ${LPG_SRC})

set_target_properties(transit_ray lpg PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_link_libraries(transit_ray PRIVATE ${GSL_LIBRARIES} ${CMAKE_DL_LIBS} omp)
target_link_libraries(lpg PRIVATE ${GSL_LIBRARIES} ${CMAKE_DL_LIBS} omp)

if(UNIX AND NOT APPLE)
  target_link_libraries(transit_ray PRIVATE m)
  target_link_libraries(lpg PRIVATE m)
endif()

# --- Optional: Custom Targets ---
set(PY_CONV_OUTPUT_DIR ${CMAKE_BINARY_DIR}/bin/conv)
add_custom_target(line_bin
    COMMAND python3 ${CMAKE_SOURCE_DIR}/conv_core/core/setup_ray.py build_ext --build-lib ${PY_CONV_OUTPUT_DIR}
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/build
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Building Python extension (ray)"
)
add_custom_target(reb_bin
    COMMAND python3 ${CMAKE_SOURCE_DIR}/conv_core/core/setup_rebin.py build_ext --build-lib ${PY_CONV_OUTPUT_DIR}
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/build
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Building Python extension (rebin)"
)
add_custom_target(all_targets DEPENDS transit_ray lpg line_bin reb_bin)

# --- Custom target to delete CMake caches ---
add_custom_target(clean_cache
  COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/CMakeFiles"
  COMMAND ${CMAKE_COMMAND} -E remove -f "${CMAKE_BINARY_DIR}/CMakeCache.txt"
  COMMAND ${CMAKE_COMMAND} -E remove -f "${CMAKE_BINARY_DIR}/cmake_install.cmake"
  COMMENT "Deleting CMake cache and CMakeFiles directory..."
)