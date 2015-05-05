file(GLOB_RECURSE CPP_FILES *.cpp)
add_library(parsa STATIC ${CPP_FILES})
