cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)
set(PROJECT "AD-PIV")
project(${PROJECT})

set(CMAKE_Fortran_FLAGS "")
set(CMAKE_Fortran_FLAGS_DEBUG "-g2 -Wall -Wno-tabs -Wno-maybe-uninitialized -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -finit-real=snan")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3")

set(common_src
	"src/array.f90"
	"src/autodiff.f90"
	"src/cluster.f90"
	"src/config.f90"
	"src/displacement.f90"
	"src/kinds.f90"
	"src/netCDF.f90"
	"src/pair.f90"
	"src/piv.f90"
	"src/settings.f90"
	"src/stats.f90"
	"src/text.f90"
	"src/time.f90"
	)

set(process_src
	"src/generator.f90"
	"src/process.f90")

find_library(NETCDF_LIBRARY netcdff "/usr/lib")
find_path(NETCDF_INCLUDE "netcdf.mod" "/usr/include")

include_directories(${NETCDF_INCLUDE})

add_custom_target(symlink
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/results" "${CMAKE_BINARY_DIR}/results"
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/figures" "${CMAKE_BINARY_DIR}/figures"
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/config" "${CMAKE_BINARY_DIR}/config"
	WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

add_executable(process ${common_src} ${process_src})
target_link_libraries(process "gomp")
target_link_libraries(process ${NETCDF_LIBRARY})
add_dependencies(process symlink)

set(FORD_EXECUTABLE "ford")
add_custom_target(ford
	COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_SOURCE_DIR}/pages" "${CMAKE_BINARY_DIR}/pages"
	COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_SOURCE_DIR}/src" "${CMAKE_BINARY_DIR}/src"
	COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/project.md" "${CMAKE_BINARY_DIR}"
	COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/README.md" "${CMAKE_BINARY_DIR}"
	COMMAND ${FORD_EXECUTABLE} "project.md"
	WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

