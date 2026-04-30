cmake_minimum_required(VERSION 3.14)
project(sharda VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Threads REQUIRED)

# ── Options ──────────────────────────────────────────────────────────────────
option(SHARDA_BUILD_TESTS "Build unit tests" ON)
option(SHARDA_FETCH_DEPS "Fetch dependencies when not available locally" ON)

# ── FetchContent dependencies ────────────────────────────────────────────────
include(FetchContent)

# spdlog
find_package(spdlog CONFIG QUIET)
if(NOT TARGET spdlog::spdlog)
    if(SHARDA_FETCH_DEPS)
        FetchContent_Declare(spdlog
            GIT_REPOSITORY https://github.com/gabime/spdlog.git
            GIT_TAG        v1.15.3)
        FetchContent_MakeAvailable(spdlog)
    else()
        message(FATAL_ERROR "spdlog not found. Install it or enable SHARDA_FETCH_DEPS")
    endif()
endif()

# HiGHS
find_library(HIGHS_LIBRARY NAMES highs HINTS ${HIGHS_DIR}/lib)
find_path(HIGHS_INCLUDE_DIR NAMES Highs.h
    HINTS ${HIGHS_DIR}/include
    PATH_SUFFIXES highs)

if(HIGHS_LIBRARY AND HIGHS_INCLUDE_DIR)
    if(NOT TARGET highs)
        add_library(highs UNKNOWN IMPORTED)
        set_target_properties(highs PROPERTIES
            IMPORTED_LOCATION "${HIGHS_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${HIGHS_INCLUDE_DIR}")
    endif()
elseif(SHARDA_FETCH_DEPS)
    FetchContent_Declare(highs
        GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
        GIT_TAG        v1.9.0)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(highs)
else()
    message(FATAL_ERROR "HiGHS not found. Install it or set -DHIGHS_DIR=<prefix>, or enable SHARDA_FETCH_DEPS")
endif()

# htslib – system installed or user provides HTSLIB_DIR
find_library(HTSLIB hts HINTS ${HTSLIB_DIR}/lib)
find_path(HTSLIB_INCLUDE htslib/sam.h HINTS ${HTSLIB_DIR}/include)
if(NOT HTSLIB OR NOT HTSLIB_INCLUDE)
    message(FATAL_ERROR "htslib not found. Install it or set -DHTSLIB_DIR=<prefix>")
endif()

# ── Library target ───────────────────────────────────────────────────────────
file(GLOB_RECURSE SHARDA_SRC src/*.cpp)
list(FILTER SHARDA_SRC EXCLUDE REGEX "src/main\\.cpp$")

add_library(sharda_lib STATIC ${SHARDA_SRC})
target_include_directories(sharda_lib PUBLIC src ${HTSLIB_INCLUDE})
target_link_libraries(sharda_lib PUBLIC spdlog::spdlog highs ${HTSLIB})

# ── Executable ───────────────────────────────────────────────────────────────
add_executable(sharda src/main.cpp)
target_link_libraries(sharda PRIVATE sharda_lib Threads::Threads)

install(TARGETS sharda RUNTIME DESTINATION bin)
install(FILES LICENSE DESTINATION share/licenses/sharda)

# ── Tests ────────────────────────────────────────────────────────────────────
if(SHARDA_BUILD_TESTS)
    FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.15.2)
    FetchContent_MakeAvailable(googletest)
    enable_testing()

    file(GLOB TEST_SRC tests/*.cpp)
    add_executable(sharda_tests ${TEST_SRC})
    target_link_libraries(sharda_tests PRIVATE sharda_lib GTest::gtest_main)
    include(GoogleTest)
    gtest_discover_tests(sharda_tests)
endif()
