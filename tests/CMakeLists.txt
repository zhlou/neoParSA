set(DIRS fly rastrigin tsp UDRST transcpp)
foreach (dir ${DIRS})
	add_subdirectory(${dir})
endforeach()
