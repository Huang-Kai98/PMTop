file(GLOB SOURCES "*.cpp")



add_library(FEM ${SOURCES})

target_include_directories(FEM PUBLIC ${CMAKE_SOURCE_DIR}/src/General  ${CMAKE_SOURCE_DIR}/src/Mesh ${CMAKE_SOURCE_DIR}/src/NumericalAlgebra)
target_link_libraries(FEM PUBLIC General Mesh NumericalAlgebra )