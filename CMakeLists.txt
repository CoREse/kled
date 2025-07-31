cmake_minimum_required(VERSION 3.15)
project(kled LANGUAGES CXX)

if(DEFINED ENV{CONDA_PREFIX})
    list(APPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
endif()

set(HTSLIB_BINARY_DIR ${CMAKE_BINARY_DIR}/htslib-build)
set(HTSLIB_INSTALL_DIR ${CMAKE_BINARY_DIR}/htslib-install)

set(HTSLIB_DEPENDENCIES "")
set(HTSLIB_DEPENDENCIES_INCLUDES "")

find_package(ZLIB 1.3 REQUIRED)
list(APPEND HTSLIB_DEPENDENCIES ${ZLIB_LIBRARIES})
list(APPEND HTSLIB_DEPENDENCIES_INCLUDES ${ZLIB_INCLUDE_DIRS})
find_package(BZip2 1.0 REQUIRED)
list(APPEND HTSLIB_DEPENDENCIES ${BZIP2_LIBRARIES})
list(APPEND HTSLIB_DEPENDENCIES_INCLUDES ${BZIP2_INCLUDE_DIRS})
find_package(LibLZMA 5.8 REQUIRED)
list(APPEND HTSLIB_DEPENDENCIES ${LIBLZMA_LIBRARIES})
list(APPEND HTSLIB_DEPENDENCIES_INCLUDES ${LIBLZMA_INCLUDE_DIRS})
find_package(CURL 8 REQUIRED)
list(APPEND HTSLIB_DEPENDENCIES ${CURL_LIBRARIES})
list(APPEND HTSLIB_DEPENDENCIES_INCLUDES ${CURL_INCLUDE_DIRS})
find_package(Threads REQUIRED)
list(APPEND HTSLIB_DEPENDENCIES ${THREADS_LIBRARIES})
list(APPEND HTSLIB_DEPENDENCIES_INCLUDES ${THREADS_INCLUDE_DIRS})
find_package(OpenSSL 3.5 REQUIRED)
list(APPEND HTSLIB_DEPENDENCIES ${OPENSSL_LIBRARIES})
list(APPEND HTSLIB_DEPENDENCIES_INCLUDES ${OPENSSL_INCLUDE_DIRS})

find_library(DEFLATE_LIB NAMES deflate)
if (NOT DEFLATE_LIB)
    message(FATAL_ERROR "libdeflate not found. Install libdeflate-dev!")
endif()

find_package(OpenMP REQUIRED)
find_package(Boost 1.84 REQUIRED COMPONENTS serialization)
find_library(GMP_LIB NAMES gmp)
if (NOT GMP_LIB)
    message(FATAL_ERROR "libgmp not found. Install gmp!")
endif()

if(UNIX)
  set(SEP ":")
else()
  set(SEP ";")
endif()

set(HTSLIB_DEPENDENCIES_ENV_L "")
foreach(p IN LISTS HTSLIB_DEPENDENCIES)
    get_filename_component(PD ${p} DIRECTORY)
endforeach()
list(JOIN HTSLIB_DEPENDENCIES_ENV_L "${SEP}" HTSLIB_DEPENDENCIES_ENV)
list(JOIN HTSLIB_DEPENDENCIES_INCLUDES "${SEP}" HTSLIB_DEPENDENCIES_INCLUDES_ENV)

include(ExternalProject)

ExternalProject_Add(
    htslib
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/htslib
    CONFIGURE_COMMAND
        sh -c "cd <SOURCE_DIR> && autoreconf -i -f && cd <BINARY_DIR> && env \
                LD_LIBRARY_PATH=\"${HTSLIB_DEPENDENCIES_ENV}:${LD_LIBRARY_PATH}\" \
                LIBRARY_PATH=\"${HTSLIB_DEPENDENCIES_ENV}:${LIBRARY_PATH}\" \
                CPATH=\"${HTSLIB_DEPENDENCIES_INCLUDES_ENV}:${CPATH}\" \
            <SOURCE_DIR>/configure --prefix=${HTSLIB_INSTALL_DIR}"
    BUILD_COMMAND $(MAKE)
    INSTALL_COMMAND $(MAKE) install
    BUILD_IN_SOURCE 0
    BUILD_BYPRODUCTS ${HTSLIB_INSTALL_DIR}/lib/libhts.a
)

add_library(hts INTERFACE)
add_dependencies(hts htslib)

target_include_directories(hts INTERFACE
    $<BUILD_INTERFACE:${HTSLIB_INSTALL_DIR}/include>
    $<INSTALL_INTERFACE:include>
    ${BZip2_INCLUDE_DIRS}
)

target_link_libraries(hts INTERFACE
    ${HTSLIB_INSTALL_DIR}/lib/libhts.a
    ${ZLIB_LIBRARIES}
    ${BZIP2_LIBRARIES}
    ${LIBLZMA_LIBRARIES}
    ${CURL_LIBRARIES}
    ${OPENSSL_LIBRARIES}
    ${DEFLATE_LIB}
    ${THREADS_LIBRARIES}
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_subdirectory(crelib)
add_subdirectory(optutils)

add_executable(kled
    kled.cpp
    input.cpp
    signature.cpp
    contig.cpp
    StatsManager.cpp
    StatsTracker.cpp
    clustering.cpp
    report.cpp
    merge.cpp
)

target_link_libraries(kled PRIVATE
    crelib_lib
    optutils_lib
    hts
    ${OpenSSL_LIBRARIES}
    m
    ${GMP_LIB}
    Boost::serialization
)

target_include_directories(kled PRIVATE
    ${HTSLIB_INSTALL_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/crelib
    ${CMAKE_CURRENT_SOURCE_DIR}/optutils
    ${BZIP2_INCLUDE_DIRS}
)

target_compile_options(kled PRIVATE -O3 -fPIE)