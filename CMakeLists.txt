cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)
set(PROJECT "main")
project(${PROJECT})

set(sources
	"src/autodiff.f90"
	"src/generator.f90"
	"src/kinds.f90"
	"src/main.f90"
	"src/netCDF.f90"
	"src/pair.f90"
	"src/piv.f90"
	"src/plplotlib.f90"
	"src/test.f90"
	"src/utilities.f90"
	)


set(CMAKE_Fortran_FLAGS "-g -Wall -Wtabs")
set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Wtabs -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -finit-real=snan")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -mtune=native -march=native")

include_directories("/usr/include")
include_directories("/usr/include/plplot")
include_directories("/usr/lib/x86_64-linux-gnu/fortran/modules/plplot")

add_executable(${PROJECT} ${sources})
target_link_libraries(${PROJECT} "plplotf95d" "plplotf95cd")
target_link_libraries(${PROJECT} "netcdff")

set(FORD_EXECUTABLE "ford")
add_custom_target(ford
	COMMAND ${FORD_EXECUTABLE} "${CMAKE_SOURCE_DIR}/project.md"
	WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
