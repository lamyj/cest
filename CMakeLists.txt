cmake_minimum_required(VERSION 3.14)

project("cest" VERSION 0.2.0)

option(BUILD_SHARED_LIBS "Build with shared libraries." ON)

set(CMAKE_INSTALL_MESSAGE LAZY)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Python COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(pybind11 REQUIRED)
find_package(SUNDIALS REQUIRED)
find_package(TBB REQUIRED)

file(
    GLOB_RECURSE cpp_files
    "${CMAKE_CURRENT_SOURCE_DIR}/cest/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/cest/*.cpp")
list(SORT cpp_files)

execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import slimp; print(slimp.__path__[0])"
    OUTPUT_VARIABLE SLIMP_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
set(compile_helper "${Python_EXECUTABLE}" "${SLIMP_PATH}/compile.py")

add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/wasabi_sampler.h"
    COMMAND ${compile_helper} compile 
        "${CMAKE_CURRENT_SOURCE_DIR}/cest/wasabi_sampler.stan"
        "${CMAKE_CURRENT_BINARY_DIR}/wasabi_sampler.h"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/cest/wasabi_sampler.stan")

pybind11_add_module(pycest ${cpp_files} "${CMAKE_CURRENT_BINARY_DIR}/wasabi_sampler.h")
target_compile_options(pycest PUBLIC -pthread -DSTAN_THREADS)
if("${TBB_VERSION}" VERSION_GREATER_EQUAL "2020")
    target_compile_options(pycest PUBLIC -DTBB_INTERFACE_NEW)
endif()
target_include_directories(
    pycest PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
    "$ENV{CMDSTAN}/stan/src" "$ENV{CMDSTAN}/stan/lib/stan_math")
target_link_libraries(
    pycest PUBLIC Eigen3::Eigen Python::NumPy SUNDIALS::generic TBB::tbb slimp)

set_target_properties(
    pycest PROPERTIES OUTPUT_NAME _cest $<$<PLATFORM_ID:Darwin>:SUFFIX .so>)

execute_process(
    COMMAND ${Python_EXECUTABLE}
        -c "import os; import sysconfig; \
            scheme = {'posix': 'posix_prefix', 'nt': 'nt'}[os.name]; \
            path = sysconfig.get_path('purelib', scheme, {'base': '.'}); \
            print(path.replace(os.path.sep, '/'))"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
    OUTPUT_STRIP_TRAILING_WHITESPACE)
    
install(DIRECTORY DESTINATION "${PYTHON_SITE_PACKAGES}")
install(TARGETS pycest DESTINATION "${PYTHON_SITE_PACKAGES}/cest")

install(
    DIRECTORY cest
    DESTINATION ${CMAKE_INSTALL_PREFIX}/${PYTHON_SITE_PACKAGES}
    FILES_MATCHING PATTERN *.py)
