cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)
set(PROJECT "AD-PIV")
project(${PROJECT})

if(CMAKE_FORTRAN_COMPILER_VERSION VERSION_GREATER 5.0)
	set(WTABS "-Wno-tabs")
else(CMAKE_FORTRAN_COMPILER_VERSION VERSION_GREATER 5.0)
	set(WTABS "-Wtabs")
endif(CMAKE_FORTRAN_COMPILER_VERSION VERSION_GREATER 5.0)

set(CMAKE_Fortran_FLAGS "")
set(CMAKE_Fortran_FLAGS_DEBUG "-g2 -Wall ${WTABS} -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -finit-real=snan")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -mtune=native -march=native")

set(common_src
	"src/autodiff.f90"
	"src/cluster.f90"
	"src/config.f90"
	"src/kinds.f90"
	"src/netCDF.f90"
	"src/pair.f90"
	"src/piv.f90"
	"src/settings.f90"
	"src/utilities.f90")

set(process_src
	"src/generator.f90"
	"src/process.f90")

set(post_src
	"src/plplotlib.f90"
	"src/post.f90")

set(figures_src
	"src/generator.f90"
	"src/plplotlib.f90"
	"src/figures.f90")

include(FindPkgConfig)
pkg_search_module(PLPLOT REQUIRED plplot-f95 plplotd-f95)
include_directories(${PLPLOT_INCLUDE_DIRS})
link_directories(${PLPLOT_LIBRARY_DIRS})

find_library(NETCDF_LIBRARY netcdff "/usr/lib")
find_path(NETCDF_INCLUDE "netcdf.mod" "/usr/include")

include_directories(${NETCDF_INCLUDE})

add_custom_target(symlink
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/results" "${CMAKE_BINARY_DIR}/results"
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/config" "${CMAKE_BINARY_DIR}/config"
	WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

add_executable(process ${common_src} ${process_src})
target_link_libraries(process "gomp")
target_link_libraries(process ${NETCDF_LIBRARY})
add_dependencies(process symlink)

add_executable(post ${common_src} ${post_src})
target_link_libraries(post "gomp")
target_link_libraries(post ${NETCDF_LIBRARY})
target_link_libraries(post ${PLPLOT_LIBRARIES})
add_dependencies(post symlink)

add_executable(figures ${common_src} ${figures_src})
target_link_libraries(figures "gomp")
target_link_libraries(figures ${NETCDF_LIBRARY})
target_link_libraries(figures ${PLPLOT_LIBRARIES})
add_dependencies(figures symlink)

set(FORD_EXECUTABLE "ford")
add_custom_target(ford
	COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_SOURCE_DIR}/pages" "${CMAKE_BINARY_DIR}/pages"
	COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_SOURCE_DIR}/src" "${CMAKE_BINARY_DIR}/src"
	COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/project.md" "${CMAKE_BINARY_DIR}"
	COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/README.md" "${CMAKE_BINARY_DIR}"
	COMMAND ${FORD_EXECUTABLE} "project.md"
	WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")