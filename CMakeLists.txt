cmake_minimum_required(VERSION 3.0)
project(Cascade VERSION 1.0)

# Include tests, which have their own CMakeLists in a subdirectory
enable_testing()
add_subdirectory(tests)

# Setup Cascade 
set(CascadeSrc
	src/cascade.c
	src/acb_ode.c
	src/examples.c
)

set(CascadeHeaders
	src/cascade.h
	src/acb_ode.h
)

add_library(cascade SHARED ${CascadeSrc})

target_link_libraries(cascade flint arb)

if(MSVC)
	target_compile_options(cascade PRIVATE /W4 /WX /O2)
else(MSVC)
	target_compile_options(cascade PRIVATE -Wall -Wextra -pedantic -O2 -g -funroll-loops)
endif(MSVC)

set_target_properties(cascade PROPERTIES PUBLIC_HEADER "${CascadeHeaders}")

install(TARGETS cascade
	LIBRARY DESTINATION lib COMPONENT devel
	PUBLIC_HEADER DESTINATION include COMPONENT devel
)
