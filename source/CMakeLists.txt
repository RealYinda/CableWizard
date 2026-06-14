cmake_minimum_required(VERSION 3.3)

project(JPSOL_example CXX C Fortran)

set(CMAKE_CXX_STANDARD 11)

set(NDIM 3)
set(CREL dbg)


set(libfind_ROOT /opt/jaumin1.14)
set(JAUMIN_ROOT  /opt/jaumin1.14)
set(JPSOL_ROOT   /opt/jaumin1.14)
list(INSERT CMAKE_MODULE_PATH 0 ${libfind_ROOT}/share/cmake/Modules/scns)
list(APPEND CMAKE_MODULE_PATH   ${JAUMIN_ROOT}/share/cmake/Modules/scns)
list(APPEND CMAKE_MODULE_PATH   ${JPSOL_ROOT}/share/cmake/Modules/scns)

set(JPSOL_BACKEND JAUMIN)
set(JAUMIN_NDIM   ${NDIM})
set(JAUMIN_CREL   ${CREL})
set(JPSOL_NDIM    ${NDIM})
set(JPSOL_CREL    ${CREL})

include(LibFindMacros)
libfind_check_external_packages(deps NAMES JPSOL ${JPSOL_BACKEND} REQUIRED)
include_directories(${deps_INCLUDE_DIRS})
add_definitions(${deps_DEFINITIONS})

include_directories(./ fem utils)
file(GLOB_RECURSE src_files ./*.C ./*.cpp ./*.c ./*.f ./*.F90
                            ./fortran/${NDIM}d/*.f
                            ./fortran/${NDIM}d/*.F90)

add_executable(main${NDIM}d ${src_files})
target_link_libraries(main${NDIM}d ${deps_LIBRARIES})
