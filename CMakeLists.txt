cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)
set(PROJECT "AD-PIV")
project(${PROJECT})

set(CMAKE_Fortran_FLAGS "-g -Wall -Wtabs -fopenmp")
set(CMAKE_Fortran_FLAGS_DEBUG "-g2 -Wall -Wtabs -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -finit-real=snan -fopenmp")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -mtune=native -march=native -fopenmp")

set(common_src
	"src/autodiff.f90"
	"src/config.f90"
	"src/generator.f90"
	"src/kinds.f90"
	"src/netCDF.f90"
	"src/pair.f90"
	"src/piv.f90"
	"src/settings.f90"
	"src/utilities.f90"
	)

set(process_src
	"src/mpiUtil.f90"
	"src/process.f90"
	)

set(post_src
	"src/plplotlib.f90"
	"src/post.f90"
	)

include_directories("/usr/include")
include_directories("/usr/include/plplot")
include_directories("/usr/lib/x86_64-linux-gnu/fortran/modules/plplot")

add_executable(process ${common_src} ${process_src})
target_link_libraries(process "gomp")
target_link_libraries(process "netcdff")

add_executable(post ${common_src} ${post_src})
target_link_libraries(post "gomp")
target_link_libraries(post "netcdff")
target_link_libraries(post "plplotf95d" "plplotf95cd")

set(FORD_EXECUTABLE "ford")
add_custom_target(ford
	COMMAND ${FORD_EXECUTABLE} "${CMAKE_SOURCE_DIR}/project.md"
	WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
