cmake_minimum_required(VERSION 3.17)
project(LammSolvF LANGUAGES Fortran)

# stdlib build options (keep build light — no quad precision, no tests)
set(CMAKE_MAXIMUM_RANK 2)
set(BUILD_TESTING OFF)
set(WITH_QP OFF)
set(WITH_XDP OFF)

# Fetch fortran-stdlib
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/stdlib/CMakeLists.txt")
    add_subdirectory(stdlib)
else()
    set(stdlib-url "https://github.com/fortran-lang/stdlib")
    message(STATUS "Retrieving stdlib from ${stdlib-url}")
    include(FetchContent)
    FetchContent_Declare(
        fortran_stdlib
        GIT_REPOSITORY ${stdlib-url}
        GIT_TAG "HEAD"
    )
    FetchContent_GetProperties(fortran_stdlib)
    if(NOT fortran_stdlib_POPULATED)
        FetchContent_Populate(fortran_stdlib)
        add_subdirectory(
            ${fortran_stdlib_SOURCE_DIR}
            ${fortran_stdlib_BINARY_DIR}
            EXCLUDE_FROM_ALL
        )
    endif()
endif()

# LAPACK and BLAS — required by the solver (dgttrf, dgttrs, ddot)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

# Shared library target
set(srcs)
add_subdirectory(src)
add_library(LammSolvF SHARED ${srcs})
target_link_libraries(LammSolvF fortran_stdlib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

# Install
if(SKBUILD)
    if(WIN32)
        set_target_properties(LammSolvF PROPERTIES LIBRARY_PREFIX "lib")
        set_target_properties(LammSolvF PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
        if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
            target_link_options(LammSolvF PRIVATE "-static")
        endif()
    endif()

    message(STATUS "Building with scikit-build")
    install(
        TARGETS LammSolvF
        LIBRARY DESTINATION lammsolvflib
        ARCHIVE DESTINATION lammsolvflib
        RUNTIME DESTINATION lammsolvflib
    )
else()
    install(
        TARGETS LammSolvF
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif()
