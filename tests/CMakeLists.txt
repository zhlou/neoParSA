set(DIRS fly rastrigin tsp)
foreach (dir ${DIRS})
	add_subdirectory(${dir})
endforeach()
