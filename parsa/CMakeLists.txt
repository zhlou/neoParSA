file(GLOB CPP_FILES *.cpp DoS/*.cpp schedule/*.cpp)
add_library(parsa STATIC ${CPP_FILES})
