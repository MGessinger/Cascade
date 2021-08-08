set(Tests
	legendre
	bessel
	hypgeom
	radius
	reduce
	solve_fuchs

	indicial_polynomial
	solve_frobenius_singleton
	solution_extend
	solution_update
	solve_frobenius
)

include_directories(${CMAKE_SOURCE_DIR}/src)

foreach (test_file IN LISTS Tests)
	add_executable(${test_file}_target ${test_file}.c)
	add_dependencies(${test_file}_target cascade)
	target_link_libraries(${test_file}_target cascade flint arb gmp)
	add_test(NAME ${test_file}_test COMMAND ${test_file}_target)
endforeach()
