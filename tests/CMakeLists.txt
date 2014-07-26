set(DIRS fly rastrigin tsp UDRST)
foreach (dir ${DIRS})
	add_subdirectory(${dir})
endforeach()
