cmake_minimum_required(VERSION 3.0)
project(Cascade VERSION 0.5.3)

set(CascadeSrc
    ../src/cascade.c
    ../src/cascade.h
    ../src/acb_ode.c
    ../src/acb_ode.h
    ../src/juliaInterface.c
    ../src/juliaInterface.h
    ../src/examples.c
    ../src/examples.h
)

add_library(cascade SHARED ${CascadeSrc})

target_link_libraries(cascade flint arb)

if(MSVC)
  target_compile_options(cascade PRIVATE /W4 /WX /O2)
else(MSVC)
  target_compile_options(cascade PRIVATE -Wall -Wextra -pedantic -O3 -s)
endif(MSVC)

set_target_properties(cascade PROPERTIES PUBLIC_HEADER "../src/cascade.h;../src/acb_ode.h;../src/examples.h")

install(TARGETS cascade
     RUNTIME DESTINATION bin COMPONENT runtime
     LIBRARY DESTINATION lib COMPONENT devel
     PUBLIC_HEADER DESTINATION include COMPONENT devel
)
