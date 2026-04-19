set(dir ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND srcs
    "${dir}/utils.f90"
    "${dir}/solver.f90"
    "${dir}/api.f90"
)
set(srcs ${srcs} PARENT_SCOPE)
