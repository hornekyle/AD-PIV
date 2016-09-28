cmake_minimum_required(VERSION 3.0)
enable_language(Fortran)

project(ADPIV
	VERSION "0.1.0")
enable_testing()

#=======================#
#= One-time Operations =#
#=======================#

find_program(EXECUTABLE_FORD "ford")

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/media")

file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
include_directories("${CMAKE_BINARY_DIR}/modules")

#=====================#
#= Set Default Flags =#
#=====================#

if(NOT HAVE_RUN)
	set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS}
		CACHE STRING "Fortran flags" FORCE)
	
	set(debug_flags "-Wall -Wno-tabs -Wno-maybe-uninitialized -Wno-unused-dummy-argument -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -finit-real=snan")
	set(CMAKE_Fortran_FLAGS_DEBUG 
		"${CMAKE_Fortran_FLAGS_DEBUG} ${debug_flags}"
		CACHE STRING "Fortran debugging flags" FORCE)
	
	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
		CACHE STRING "Fortran release flags" FORCE)
endif(NOT HAVE_RUN)

#===============================#
#= Find Packages and Libraries =#
#===============================#

# find_file(CommonModules-cmake "CommonModules.cmake"
# 	HINTS "/usr/lib/CommonModules")
# include_directories("/usr/include/CommonModules")

find_file(CommonModules-cmake "CommonModules.cmake"
	HINTS "/home/kyle/Development/Fortran/CommonModules/build")
include(${CommonModules-cmake})
include_directories("/home/kyle/Development/Fortran/CommonModules/build/modules")

#=========================#
#= Create Symlink Target =#
#=========================#

add_custom_target(symlinks
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/results" "${CMAKE_BINARY_DIR}/results"
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/figures" "${CMAKE_BINARY_DIR}/figures"
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/config" "${CMAKE_BINARY_DIR}/config"
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/src"   "${CMAKE_BINARY_DIR}/src"
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/pages" "${CMAKE_BINARY_DIR}/pages"
	COMMAND ${CMAKE_COMMAND} -E create_symlink "${CMAKE_SOURCE_DIR}/project.md" "${CMAKE_BINARY_DIR}/project.md"
	WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

#=========================#
#= Create Library Target =#
#=========================#

# Source Files
set(process_src
	"src/cluster.f90"
	"src/displacement.f90"
	"src/generator.f90"
	"src/pair.f90"
	"src/piv.f90"
	"src/process.f90"
	"src/settings.f90"
	)

# Add target
add_executable(process ${process_src})
target_link_libraries(process CommonModules)
add_dependencies(process symlinks)

#======================#
#= Create FORD Target =#
#======================#

add_custom_command(OUTPUT "docs/index.html"
	COMMAND basic
	COMMAND examples
	COMMAND ${EXECUTABLE_FORD} "${CMAKE_BINARY_DIR}/project.md"
	DEPENDS logo basic examples symlinks
	WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
add_custom_target(doc DEPENDS "docs/index.html")

set(HAVE_RUN "1" CACHE STRING "Indicate if first run of CMake" FORCE)
