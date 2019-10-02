cmake_minimum_required(VERSION 3.0)
project(Cascade VERSION 1.0)

set(CascadeSrc
    src/cascade.c
    src/cascade.h
    src/acb_ode.c
    src/acb_ode.h
    src/examples.c
)

add_library(cascade SHARED ${CascadeSrc})

target_link_libraries(cascade flint arb)

if(MSVC)
  target_compile_options(cascade PRIVATE /W4 /WX /O2)
else(MSVC)
  target_compile_options(cascade PRIVATE -Wall -Wextra -pedantic -O3)
endif(MSVC)

set_target_properties(cascade PROPERTIES PUBLIC_HEADER "src/cascade.h;src/acb_ode.h")

install(TARGETS cascade
     LIBRARY DESTINATION lib COMPONENT devel
     PUBLIC_HEADER DESTINATION include COMPONENT devel
)
