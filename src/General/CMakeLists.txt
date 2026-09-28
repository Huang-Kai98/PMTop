file(GLOB SOURCES "*.cpp")

add_library(General ${SOURCES})
target_include_directories(General PRIVATE ${CMAKE_SOURCE_DIR}/src/NumericalAlgebra)
target_link_libraries(General PRIVATE NumericalAlgebra) 
