cmake_minimum_required(VERSION 3.10.2)
project(ee382n_bitcoin)

enable_language(C ASM)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
find_package(CURL REQUIRED)

add_subdirectory(compat)

add_compile_options(-fno-strict-aliasing)

set(minerd_SOURCES
        elist.h
        miner.h
        compat.h
        cpu-miner.c
        util.c
        sha2.c
        scrypt.c
        cpuminer-config.h
        # Architecture-specific assembly hash functions
        sha2-arm.S
        sha2-ppc.S
        sha2-x64.S
        sha2-x86.S
        scrypt-arm.S
        scrypt-ppc.S
        scrypt-x64.S
        scrypt-x86.S
        # ECE 382N: FPGA driver files
        fpga_sha.c
        fpga_sha.h
)

add_executable(minerd ${minerd_SOURCES})
include_directories(${CURL_INCLUDE_DIR})
target_link_libraries(minerd
#        CURL::libcurl
        ${CURL_LIBRARIES}
        Threads::Threads
        jansson
)
