cmake_minimum_required(VERSION 3.14)
project(multiply_gates)

# Example usages:
# 1. Run multiply_gates tests (requires Catch2v3):
#   mkdir build && cd build && cmake -DBUILD_TESTS=ON .. && cmake --build . \
#     && ./tests/test_multiply_gates
# 2. Build multiply_gates shared library:
#   mkdir build && cd build \
#     && cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Debug .. && cmake --build .
# 3. Release build multiply_gates shared library:
#   mkdir build && cd build \
#     && cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release .. \
#     && cmake --build .
# 4. Cross-platform release build multiply_gates shared library on Linux targeting Windows:
#   mkdir build && cd build \
#     && cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release \
#       -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchains/x86_64-w64-mingw32.cmake .. \
#     && cmake --build .
# 5. You may want to copy that library after compilation:
#   # assuming you are in build directory and ran #2 or #3:
#   cp src/multiply_gates/cpp_src/libmultiply_gates.so ../src/multiply_gates/lib/
# 6. Cleanup
#   cd .. # if you are in build directory
#   rm -r build

option(BUILD_TESTS "Build tests" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

if(BUILD_SHARED_LIBS)
  add_subdirectory(src/genftqc/cpp_src)
  # If ${SKBUILD_PLATLIB_DIR} is defined, we need to install the shared library
  # (compiled in src/genftqc/cpp_src) to that directory:
  if(DEFINED SKBUILD_PLATLIB_DIR)
    install(TARGETS genftqc
      LIBRARY DESTINATION ${SKBUILD_PLATLIB_DIR}/genftqc/lib
      RUNTIME DESTINATION ${SKBUILD_PLATLIB_DIR}/genftqc/lib
    )
  endif()
endif()

if(BUILD_TESTS)
  add_subdirectory(tests)
endif()
