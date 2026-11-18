cmake_minimum_required(VERSION 3.13)

file(STRINGS include/codec2/version.h _ver REGEX "^#define C2FXP_VERSION ")
string(REGEX REPLACE ".*\"(.*)\".*" "\\1" C2FXP_VERSION "${_ver}")
project(codec2-fixed-point VERSION ${C2FXP_VERSION} LANGUAGES C)

include(GNUInstallDirs)

option(BUILD_SHARED_LIBS "Build libcodec2 as a shared library" ON)
option(C2FXP_BUILD_TOOLS "Build the c2fxp command-line tool" ON)
option(C2FXP_BUILD_TESTS "Build the tests (ctest)" ON)
option(C2FXP_SANITIZE "Build the tests with UBSan, ASan and the FXP_CHECK bound asserts" OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CODEC_SRCS
    src/analysis.c src/codec2_fxp.c src/delta_lsp_cb.c src/fxp_analysis_tables.c src/fxp_fft.c
    src/fxp_fft_tables.c src/fxp_lpc_tables.c src/fxp_lpcspec.c src/fxp_math.c src/fxp_math_tables.c
    src/fxp_nlp_tables.c src/fxp_rdft.c src/interp.c src/lpc.c src/nlp.c src/quantise.c src/synthesis.c)

# libcodec2: the integer codec plus the upstream-compatible API, under upstream's library name
add_library(codec2 ${CODEC_SRCS} src/codec2_compat.c)
target_include_directories(codec2
    PRIVATE src
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include> $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/codec2>
    INTERFACE $<INSTALL_INTERFACE:include/codec2> $<INSTALL_INTERFACE:include>)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(codec2 PRIVATE -Wall -Wextra)
endif()
if(BUILD_SHARED_LIBS)
    target_compile_definitions(codec2 PRIVATE C2FXP_BUILDING_SHARED)
    # export only the public API; the codec's internal functions stay private to the library
    set_target_properties(codec2 PROPERTIES C_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
endif()
# same soname as upstream Codec2 1.2 (libcodec2.so.1.2), so it can stand in for it at run time
set_target_properties(codec2 PROPERTIES SOVERSION 1.2
    PUBLIC_HEADER "include/codec2/codec2.h;include/codec2/codec2_fxp.h;include/codec2/version.h")

install(TARGETS codec2 EXPORT codec2-config
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/codec2)
# as upstream: find_package(codec2) gives the imported target "codec2"
install(EXPORT codec2-config DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/codec2)
export(TARGETS codec2 FILE ${PROJECT_BINARY_DIR}/codec2.cmake)
configure_file(cmake/codec2.pc.in ${PROJECT_BINARY_DIR}/codec2.pc @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/codec2.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

if(C2FXP_BUILD_TOOLS)
    add_executable(c2fxp tools/c2fxp.c)
    target_link_libraries(c2fxp PRIVATE codec2)
    if(APPLE)
        set_target_properties(c2fxp PROPERTIES INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
    elseif(UNIX)
        set_target_properties(c2fxp PROPERTIES INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
    endif()
    install(TARGETS c2fxp RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

if(C2FXP_BUILD_TESTS)
    include(CTest)
    find_library(MATH_LIB m)
    set(TEST_FLAGS)
    if(C2FXP_SANITIZE)
        set(TEST_FLAGS -fsanitize=undefined,address -fno-sanitize-recover=all)
    endif()

    # public API and bit-exact regression, against the library as installed users get it
    add_executable(test_api tests/test_api.c)
    target_link_libraries(test_api PRIVATE codec2 ${MATH_LIB})
    add_test(NAME api COMMAND test_api)

    # kernel unit tests: built from the codec sources with the FXP_CHECK bound asserts on
    function(c2fxp_unit name)
        add_executable(${name} tests/${name}.c ${ARGN})
        target_include_directories(${name} PRIVATE src include tests/ref)
        target_compile_definitions(${name} PRIVATE FXP_CHECK)
        target_compile_options(${name} PRIVATE ${TEST_FLAGS})
        target_link_options(${name} PRIVATE ${TEST_FLAGS})
        target_link_libraries(${name} PRIVATE ${MATH_LIB})
        add_test(NAME ${name} COMMAND ${name})
    endfunction()
    set(KISS tests/ref/kiss_fft.c tests/ref/kiss_fftr.c)
    c2fxp_unit(test_fxp_fft src/fxp_fft.c src/fxp_rdft.c src/fxp_fft_tables.c ${KISS})
    c2fxp_unit(test_fxp_lpcspec src/fxp_lpcspec.c src/fxp_math.c src/fxp_math_tables.c src/fxp_rdft.c
        src/fxp_fft_tables.c src/delta_lsp_cb.c)
    c2fxp_unit(test_decode_params src/quantise.c src/fxp_lpc_tables.c src/interp.c src/delta_lsp_cb.c)
    c2fxp_unit(test_analysis src/nlp.c src/fxp_nlp_tables.c src/fxp_analysis_tables.c src/fxp_fft.c
        src/fxp_fft_tables.c src/fxp_math.c src/fxp_math_tables.c ${KISS})
    c2fxp_unit(test_nlp src/fxp_nlp_tables.c src/fxp_fft.c src/fxp_fft_tables.c src/fxp_math.c
        src/fxp_math_tables.c ${KISS})
    c2fxp_unit(test_lpc src/analysis.c src/codec2_fxp.c src/delta_lsp_cb.c src/fxp_analysis_tables.c
        src/fxp_fft.c src/fxp_fft_tables.c src/fxp_lpc_tables.c src/fxp_lpcspec.c src/fxp_math.c
        src/fxp_math_tables.c src/fxp_nlp_tables.c src/fxp_rdft.c src/interp.c src/nlp.c src/quantise.c
        src/synthesis.c ${KISS})
endif()
