cmake_minimum_required(VERSION 3.0)
project(Cascade VERSION 0.5.3)

set(CascadeSrc
    src/cascade.c
    src/acb_ode.c
    src/juliaInterface.c
)

add_library(cascade SHARED ${CascadeSrc})

target_link_libraries(cascade flint arb gmp)

if(MSVC)
  target_compile_options(cascade PRIVATE /W4 /WX /O2)
else(MSVC)
  target_compile_options(cascade PRIVATE -Wall -Wextra -pedantic -O3 -s)
endif(MSVC)
