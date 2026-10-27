#!/usr/bin/env bash
#
set -euo pipefail
VCPKG="${VCPKG_ROOT:?set VCPKG_ROOT to your vcpkg checkout}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"       # repo root
OVL="$ROOT/native/ios-build"
BUILD="$ROOT/native/build_iossim"
DEPS="${IOS_DEPS_DIR:-$BUILD/vcpkg_installed}"      # where the ios-sim deps install
MANIFEST="${IOS_MANIFEST_DIR:?set IOS_MANIFEST_DIR to a manifest dir requiring operon (+ foolnotion registry + fork overlays)}"

echo ">> [1/2] cross-building operon dep stack for arm64-ios-simulator"
"$VCPKG/vcpkg" install \
  --x-manifest-root="$MANIFEST" \
  --x-install-root="$DEPS" \
  --overlay-triplets="$OVL/overlay-triplets" \
  --overlay-ports="$OVL/overlay-ports" \
  --triplet arm64-ios-simulator \
  --host-triplet arm64-osx \
  --no-print-usage

echo ">> [2/2] configuring + building sketch_bridge"
cmake -S "$ROOT/native" -B "$BUILD" \
  -DCMAKE_TOOLCHAIN_FILE="$VCPKG/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_TARGET_TRIPLET=arm64-ios-simulator \
  -DVCPKG_OVERLAY_TRIPLETS="$OVL/overlay-triplets" \
  -DVCPKG_OVERLAY_PORTS="$OVL/overlay-ports" \
  -DVCPKG_INSTALLED_DIR="$DEPS" \
  -DSKETCH_TARGET_PLATFORM=ios-simulator \
  -DCMAKE_SYSTEM_NAME=iOS \
  -DCMAKE_OSX_SYSROOT=iphonesimulator \
  -DCMAKE_OSX_ARCHITECTURES=arm64 \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_DEPLOYMENT_TARGET="${IOS_MIN:-17.0}" \
  -DSKETCH_LTO_NATIVE=OFF
# compiled at runtime from source, so no .metallib goes in the bundle). test_bridge_gpu_pointseg is a
# C-ABI harness; run it with `xcrun simctl spawn booted <dir>/test_bridge_gpu_pointseg counters|fingerprint|solve`
# from an INTERNAL-disk dir holding it + both dylibs, with any /Volumes/... LC_RPATH deleted
# (dyld_sim blocks opening files on the external drive) and @loader_path added.
cmake --build "$BUILD" --config Release --target sketch_bridge test_bridge_gpu_pointseg -j

echo ">> done. verify:"
echo "   lipo -info $BUILD/libsketch_bridge.dylib"
echo "   otool -l $BUILD/libsketch_bridge.dylib | grep -A4 LC_BUILD_VERSION   # expect platform 7 (iOS Simulator)"
# For a physical device: use triplet arm64-ios + CMAKE_OSX_SYSROOT=iphoneos + a matching overlay
# triplet (arm64-ios); the same asmjit overlay-port applies. Device also needs code-signing (user identity).

# --- Runtime dep embedding (infix-parser) — libsketch_bridge dlopen needs its sibling dylib ---
# After building the ios-sim engine, embed libinfix-parser for the app bundle:
#   cp <deps>/arm64-ios-simulator/lib/libinfix-parser.0.dylib ios/GraphArt/libinfix-parser.0.dylib
#   install_name_tool -id "@rpath/libinfix-parser.0.dylib" ios/GraphArt/libinfix-parser.0.dylib
#   install_name_tool -delete_rpath <stale-abs-rpath> -add_rpath @loader_path ios/GraphArt/libsketch_bridge.dylib
# The Xcode project references both dylibs in Copy Bundle Resources (project.pbxproj); @loader_path lets
# sketch_bridge resolve the sibling infix-parser inside the flat iOS app bundle.
