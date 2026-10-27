# Our fork of Operon with joint t-optimization support
set(SOURCE_PATH "$ENV{OPERON_SRC}")

set(VCPKG_BUILD_TYPE release)

string(COMPARE EQUAL "${VCPKG_LIBRARY_LINKAGE}" "dynamic" OPERON_SHARED_LIBS)

vcpkg_cmake_configure(
  SOURCE_PATH "${SOURCE_PATH}"
  OPTIONS
        -DBUILD_TESTING=OFF
        -DBUILD_EXAMPLES=OFF
        -DBUILD_SHARED_LIBS=${OPERON_SHARED_LIBS}
        -DBUILD_CLI_PROGRAMS=OFF
        -DUSE_SINGLE_PRECISION=ON
        -DUSE_ASMJIT=OFF
  MAYBE_UNUSED_VARIABLES
        BUILD_TESTING
        USE_ASMJIT
)
vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME operon DO_NOT_DELETE_PARENT_CONFIG_PATH)

# asmjit is not usable on iOS; drop it from the installed config.
set(_operon_config "${CURRENT_PACKAGES_DIR}/share/operon/operonConfig.cmake")
if(EXISTS "${_operon_config}")
  file(READ "${_operon_config}" _operon_config_contents)
  string(REPLACE "find_dependency(asmjit CONFIG REQUIRED)\n" "" _operon_config_contents "${_operon_config_contents}")
  file(WRITE "${_operon_config}" "${_operon_config_contents}")
endif()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include"
                    "${CURRENT_PACKAGES_DIR}/debug/share")

file(
  INSTALL "${SOURCE_PATH}/LICENSE"
  DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}"
  RENAME copyright)
