cmake_minimum_required(VERSION 3.24)
project(graphiti-engine CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Compiler launcher: ccache when found.
find_program(SKETCH_CCACHE_PROGRAM ccache)
option(SKETCH_CCACHE "Use ccache as the compiler launcher, if found" ON)
if(SKETCH_CCACHE AND SKETCH_CCACHE_PROGRAM)
    set(CMAKE_C_COMPILER_LAUNCHER "${SKETCH_CCACHE_PROGRAM}")
    set(CMAKE_CXX_COMPILER_LAUNCHER "${SKETCH_CCACHE_PROGRAM}")
    set(CMAKE_OBJCXX_COMPILER_LAUNCHER "${SKETCH_CCACHE_PROGRAM}")
    message(STATUS "sketch: ccache enabled (${SKETCH_CCACHE_PROGRAM})")
elseif(SKETCH_CCACHE)
    message(STATUS "sketch: ccache requested but not found; building without it")
endif()

# Target platform: macos, ios-device, or ios-simulator (inferred from SDK when empty).
set(SKETCH_TARGET_PLATFORM "" CACHE STRING
    "Bridge target: macos, ios-device, or ios-simulator (empty = infer from SDK)")
if(SKETCH_TARGET_PLATFORM STREQUAL "")
    string(TOLOWER "${CMAKE_OSX_SYSROOT}" _sketch_sysroot_lower)
    if(_sketch_sysroot_lower MATCHES "iphoneos")
        set(SKETCH_TARGET_PLATFORM "ios-device")
    elseif(_sketch_sysroot_lower MATCHES "iphonesimulator")
        set(SKETCH_TARGET_PLATFORM "ios-simulator")
    else()
        set(SKETCH_TARGET_PLATFORM "macos")
    endif()
endif()
set_property(CACHE SKETCH_TARGET_PLATFORM PROPERTY STRINGS macos ios-device ios-simulator)
if(NOT SKETCH_TARGET_PLATFORM STREQUAL "macos" AND
   NOT SKETCH_TARGET_PLATFORM STREQUAL "ios-device" AND
   NOT SKETCH_TARGET_PLATFORM STREQUAL "ios-simulator")
    message(FATAL_ERROR "SKETCH_TARGET_PLATFORM must be macos, ios-device, or ios-simulator")
endif()
set(_sketch_is_ios FALSE)
if(SKETCH_TARGET_PLATFORM MATCHES "^ios-")
    set(_sketch_is_ios TRUE)
endif()
message(STATUS "sketch: target platform = ${SKETCH_TARGET_PLATFORM}")

# LTO and host/iOS CPU tuning. No fast-math.
option(SKETCH_LTO_NATIVE "LTO + -mcpu=native for the native engine (no fast-math)" ON)
if(SKETCH_LTO_NATIVE)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT _sketch_ipo_ok OUTPUT _sketch_ipo_msg)
    if(_sketch_ipo_ok)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
        message(STATUS "sketch: LTO (IPO) enabled")
    else()
        message(STATUS "sketch: LTO not supported here: ${_sketch_ipo_msg}")
    endif()
    if(NOT _sketch_is_ios)
        include(CheckCXXCompilerFlag)
        check_cxx_compiler_flag("-mcpu=native" _sketch_mcpu_native)
        if(_sketch_mcpu_native)
            add_compile_options(-mcpu=native)
            message(STATUS "sketch: -mcpu=native enabled")
        endif()
    else()
        # iOS deployment target 17.0 minimum device is A12.
        include(CheckCXXCompilerFlag)
        check_cxx_compiler_flag("-mcpu=apple-a12" _sketch_mcpu_ios)
        if(_sketch_mcpu_ios)
            add_compile_options(-mcpu=apple-a12)
            message(STATUS "sketch: -mcpu=apple-a12 enabled for iOS (derived from IPHONEOS_DEPLOYMENT_TARGET 17.0)")
        else()
            message(STATUS "sketch: -mcpu=apple-a12 not supported by this compiler; omitting (iOS gets generic arm64 codegen)")
        endif()
    endif()
endif()

# PGO: 'generate' or 'use' with a merged .profdata. Empty means no PGO.
set(SKETCH_PGO "" CACHE STRING "PGO mode: '', 'generate', or 'use'")
set(SKETCH_PGO_PROFILE "" CACHE FILEPATH "merged .profdata for SKETCH_PGO=use")
if(SKETCH_PGO STREQUAL "generate")
    add_compile_options(-fprofile-instr-generate)
    add_link_options(-fprofile-instr-generate)
    message(STATUS "sketch: PGO instrumentation build (SKETCH_PGO=generate)")
elseif(SKETCH_PGO STREQUAL "use")
    if(SKETCH_PGO_PROFILE STREQUAL "")
        message(FATAL_ERROR "SKETCH_PGO=use requires -DSKETCH_PGO_PROFILE=<merged .profdata>")
    endif()
    add_compile_options(-fprofile-instr-use=${SKETCH_PGO_PROFILE}
                         -Wno-profile-instr-out-of-date -Wno-profile-instr-unprofiled)
    add_link_options(-fprofile-instr-use=${SKETCH_PGO_PROFILE})
    message(STATUS "sketch: PGO-optimized build using ${SKETCH_PGO_PROFILE}")
elseif(NOT SKETCH_PGO STREQUAL "")
    message(FATAL_ERROR "SKETCH_PGO must be '', 'generate', or 'use' (got '${SKETCH_PGO}')")
endif()

# mimalloc: optional drop-in allocator, never fails the build.
option(SKETCH_MIMALLOC "link mimalloc (drop-in allocator)" ON)
set(SKETCH_MIMALLOC_TARGET "")
if(SKETCH_MIMALLOC AND NOT _sketch_is_ios)
    find_package(mimalloc CONFIG QUIET)
    if(mimalloc_FOUND)
        set(SKETCH_MIMALLOC_TARGET "mimalloc")
        message(STATUS "sketch: mimalloc found via find_package (CONFIG)")
    else()
        execute_process(
            COMMAND brew --prefix mimalloc
            OUTPUT_VARIABLE _sketch_mimalloc_prefix
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
            RESULT_VARIABLE _sketch_brew_result
        )
        if(_sketch_brew_result EQUAL 0 AND EXISTS "${_sketch_mimalloc_prefix}/lib")
            file(GLOB _sketch_mimalloc_dylib "${_sketch_mimalloc_prefix}/lib/libmimalloc*.dylib")
            if(_sketch_mimalloc_dylib)
                list(GET _sketch_mimalloc_dylib 0 SKETCH_MIMALLOC_TARGET)
                message(STATUS "sketch: mimalloc found via Homebrew: ${SKETCH_MIMALLOC_TARGET}")
            endif()
        elseif(EXISTS "/opt/homebrew/opt/mimalloc/lib")
            file(GLOB _sketch_mimalloc_dylib "/opt/homebrew/opt/mimalloc/lib/libmimalloc*.dylib")
            if(_sketch_mimalloc_dylib)
                list(GET _sketch_mimalloc_dylib 0 SKETCH_MIMALLOC_TARGET)
                message(STATUS "sketch: mimalloc found at /opt/homebrew/opt/mimalloc: ${SKETCH_MIMALLOC_TARGET}")
            endif()
        endif()
        if(SKETCH_MIMALLOC_TARGET STREQUAL "")
            message(STATUS "sketch: mimalloc not found (find_package + Homebrew probes both empty) — skipping, build continues without it")
        endif()
    endif()
elseif(SKETCH_MIMALLOC)
    message(STATUS "sketch: mimalloc disabled for iOS target (no target-correct allocator configured)")
endif()

# Package roots: all vcpkg packages (operon and dependencies) under VCPKG_PKG_ROOT.
if(NOT DEFINED VCPKG_PKG_ROOT)
    set(VCPKG_PKG_ROOT "$ENV{VCPKG_ROOT}/packages")
endif()
if(NOT DEFINED SKETCH_VCPKG_TRIPLET)
    if(SKETCH_TARGET_PLATFORM STREQUAL "ios-device")
        set(SKETCH_VCPKG_TRIPLET "arm64-ios")
    elseif(SKETCH_TARGET_PLATFORM STREQUAL "ios-simulator")
        set(SKETCH_VCPKG_TRIPLET "arm64-ios-simulator")
    else()
        set(SKETCH_VCPKG_TRIPLET "arm64-osx")
    endif()
endif()
set(SKETCH_VCPKG_TRIPLET "${SKETCH_VCPKG_TRIPLET}" CACHE STRING
    "vcpkg package triplet matching SKETCH_TARGET_PLATFORM")
file(GLOB VCPKG_PKG_DIRS "${VCPKG_PKG_ROOT}/*_${SKETCH_VCPKG_TRIPLET}")
if(VCPKG_PKG_DIRS STREQUAL "")
    message(FATAL_ERROR
        "No packages for '${SKETCH_VCPKG_TRIPLET}' under ${VCPKG_PKG_ROOT}. "
        "A ${SKETCH_TARGET_PLATFORM} bridge requires target-correct Operon and dependencies; "
        "the macOS arm64 packages must not be reused. Configure with "
        "-DSKETCH_VCPKG_TRIPLET=<target triplet> after installing those packages.")
endif()
foreach(PKG_DIR ${VCPKG_PKG_DIRS})
    list(APPEND CMAKE_PREFIX_PATH "${PKG_DIR}")
endforeach()

find_package(operon CONFIG REQUIRED)

set(VCPKG_ALL_INCLUDES "")
foreach(PKG_DIR ${VCPKG_PKG_DIRS})
    if(EXISTS "${PKG_DIR}/include")
        list(APPEND VCPKG_ALL_INCLUDES "${PKG_DIR}/include")
    endif()
endforeach()

# Header-only engine: include dirs, operon, mimalloc, Accelerate.
add_library(sketch-engine INTERFACE)
target_include_directories(sketch-engine INTERFACE include)
target_include_directories(sketch-engine SYSTEM INTERFACE ${VCPKG_ALL_INCLUDES})
target_link_libraries(sketch-engine INTERFACE operon::operon)
if(NOT SKETCH_MIMALLOC_TARGET STREQUAL "")
    target_link_libraries(sketch-engine INTERFACE "${SKETCH_MIMALLOC_TARGET}")
endif()
if(APPLE)
    find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
    add_compile_definitions(SKETCH_HAS_ACCELERATE)
    target_link_libraries(sketch-engine INTERFACE "${ACCELERATE_FRAMEWORK}")
    message(STATUS "sketch: Apple Accelerate framework linked for hardware vectorization")
endif()

# Python binding (optional, requires nanobind).
find_package(Python COMPONENTS Interpreter Development.Module QUIET)
if(Python_FOUND)
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -c "import nanobind; print(nanobind.cmake_dir())"
        OUTPUT_VARIABLE NB_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE NB_RESULT
    )
    if(NB_RESULT EQUAL 0)
        list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
        find_package(nanobind CONFIG REQUIRED)

        nanobind_add_module(sketch_native src/python_binding.cpp)
        target_link_libraries(sketch_native PRIVATE sketch-engine)
        message(STATUS "Python binding enabled (sketch_native)")
    else()
        message(STATUS "nanobind not found — Python binding disabled")
    endif()
endif()

# Metal runtime static libraries (macOS).
if(APPLE)
    enable_language(OBJCXX)

    add_library(sketch-fp16-gram-runtime-metal STATIC src/fp16_gram_runtime_metal.mm)
    target_include_directories(sketch-fp16-gram-runtime-metal PUBLIC include)
    target_compile_options(sketch-fp16-gram-runtime-metal PRIVATE -fobjc-arc)
    target_link_libraries(sketch-fp16-gram-runtime-metal PRIVATE "-framework Metal" "-framework Foundation")

    add_library(sketch-bracket-newton-metal STATIC src/bracket_newton_runtime_metal.mm)
    target_include_directories(sketch-bracket-newton-metal PUBLIC include)
    target_link_libraries(sketch-bracket-newton-metal PRIVATE "-framework Metal" "-framework Foundation")
endif()

# C bridge (product library).
add_library(sketch_bridge SHARED src/sketch_c_bridge.cpp)
target_include_directories(sketch_bridge PUBLIC include)
target_include_directories(sketch_bridge SYSTEM PUBLIC ${VCPKG_ALL_INCLUDES})
target_link_libraries(sketch_bridge PRIVATE operon::operon sketch-engine)
if(APPLE)
    target_link_libraries(sketch_bridge PRIVATE "${ACCELERATE_FRAMEWORK}")
endif()

# Bridge built as Objective-C++ (ARC) with Metal point-seg scoring. The MSL kernel is generated
# from shaders/pointseg.metal at configure time. Set SKETCH_BRIDGE_METAL=OFF for plain C++.
option(SKETCH_BRIDGE_METAL "Build sketch_bridge as Objective-C++ with the Metal point-seg path" ON)
if(APPLE AND SKETCH_BRIDGE_METAL)
    enable_language(OBJCXX)

    # Embed shaders/pointseg.metal as a C++ raw-string header (single source of truth).
    set(SKETCH_POINTSEG_METAL_SRC "${CMAKE_CURRENT_SOURCE_DIR}/shaders/pointseg.metal")
    set(SKETCH_POINTSEG_METAL_GEN_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated")
    set(SKETCH_POINTSEG_METAL_GEN_HEADER "${SKETCH_POINTSEG_METAL_GEN_DIR}/sketch/pointseg_metal_source.hpp")
    file(MAKE_DIRECTORY "${SKETCH_POINTSEG_METAL_GEN_DIR}/sketch")
    execute_process(
        COMMAND ${CMAKE_COMMAND}
                -DIN_FILE=${SKETCH_POINTSEG_METAL_SRC}
                -DOUT_FILE=${SKETCH_POINTSEG_METAL_GEN_HEADER}
                -DVAR_NAME=kShaderGenerated
                -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/EmbedTextFile.cmake"
        RESULT_VARIABLE _sketch_embed_result
    )
    if(NOT _sketch_embed_result EQUAL 0)
        message(FATAL_ERROR "sketch: failed to generate ${SKETCH_POINTSEG_METAL_GEN_HEADER} from ${SKETCH_POINTSEG_METAL_SRC}")
    endif()
    include_directories("${SKETCH_POINTSEG_METAL_GEN_DIR}")

    # Precompiled pointseg.metallib (no fast-math, matching the runtime MTLMathModeSafe).
    set(SKETCH_POINTSEG_METALLIB "")
    find_program(SKETCH_XCRUN_EXE xcrun)
    if(SKETCH_XCRUN_EXE)
        set(_air "${CMAKE_CURRENT_BINARY_DIR}/pointseg.air")
        set(_lib "${CMAKE_CURRENT_BINARY_DIR}/pointseg.metallib")
        execute_process(
            COMMAND ${SKETCH_XCRUN_EXE} -sdk macosx metal -fno-fast-math -c "${SKETCH_POINTSEG_METAL_SRC}" -o "${_air}"
            RESULT_VARIABLE _sketch_metal_cc_result
            ERROR_VARIABLE _sketch_metal_cc_err
        )
        if(_sketch_metal_cc_result EQUAL 0)
            execute_process(
                COMMAND ${SKETCH_XCRUN_EXE} -sdk macosx metallib "${_air}" -o "${_lib}"
                RESULT_VARIABLE _sketch_metallib_result
                ERROR_VARIABLE _sketch_metallib_err
            )
            if(_sketch_metallib_result EQUAL 0)
                set(SKETCH_POINTSEG_METALLIB "${_lib}")
                message(STATUS "sketch: precompiled pointseg.metallib at ${_lib} (SKETCH_METALLIB=1 to use it)")
            else()
                message(WARNING "sketch: xcrun metallib failed: ${_sketch_metallib_err}")
            endif()
        else()
            message(WARNING "sketch: xcrun metal -c failed: ${_sketch_metal_cc_err}")
        endif()
    else()
        message(WARNING "sketch: xcrun not found -- pointseg.metallib not built, SKETCH_METALLIB=1 will fall back to newLibraryWithSource at runtime")
    endif()
    if(SKETCH_POINTSEG_METALLIB)
        add_compile_definitions(SKETCH_POINTSEG_METALLIB_PATH="${SKETCH_POINTSEG_METALLIB}")
    endif()
    set_source_files_properties(src/sketch_c_bridge.cpp TARGET_DIRECTORY sketch_bridge PROPERTIES LANGUAGE OBJCXX)
    set_target_properties(sketch_bridge PROPERTIES OBJCXX_STANDARD 23 OBJCXX_STANDARD_REQUIRED ON)
    target_compile_options(sketch_bridge PRIVATE $<$<COMPILE_LANGUAGE:OBJCXX>:-fobjc-arc>)
    target_link_libraries(sketch_bridge PRIVATE "-framework Metal" "-framework Foundation")
    message(STATUS "sketch: sketch_bridge built as Objective-C++ with Metal point-seg scoring")
endif()
