set(DIRS fly rastrigin tsp UDRST transcpp schwefel)
foreach (dir ${DIRS})
	add_subdirectory(${dir})
endforeach()
