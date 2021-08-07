cmake_minimum_required(VERSION 3.0)
project(Cascade VERSION 1.0 DESCRIPTION "Cascade - the C Library for Approximative Solutions to Complex Arbitrary precision Differential Equations")

set(DEPS gmp flint arb)

foreach (LIB ${DEPS})
	find_library(${LIB}_LIBRARY NAMES ${LIB})
	if (NOT ${LIB}_LIBRARY)
		message(FATAL_ERROR "~~ Library ${LIB} was not found.")
	endif()
	add_library(${LIB} UNKNOWN IMPORTED)
	set_property(TARGET ${LIB} PROPERTY IMPORTED_LOCATION ${${LIB}_LIBRARY})
	message("-- Library ${LIB} was found in ${${LIB}_LIBRARY}")
endforeach ()

# Include tests, which have their own CMakeLists in a subdirectory
enable_testing()
add_subdirectory(tests)

# Setup Cascade 
set(CascadeSrc
	src/acb_ode.c
	src/acb_ode_solution.c
	src/examples.c
	src/solver_fuchs.c
	src/solver_frobenius.c
)

set(CascadeHeaders
	src/acb_ode.h
	src/cascade.h
)

add_library(cascade SHARED ${CascadeSrc})
target_link_libraries(cascade ${DEPS})
target_compile_options(cascade PRIVATE -Wall -Wextra -pedantic)

set_target_properties(cascade PROPERTIES PUBLIC_HEADER "${CascadeHeaders}")
set_target_properties(cascade PROPERTIES VERSION ${PROJECT_VERSION})

install(TARGETS cascade LIBRARY PUBLIC_HEADER)
