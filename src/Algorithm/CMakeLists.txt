file(GLOB SOURCES "*.cpp")

add_library(Algorithm ${SOURCES})
target_include_directories(Algorithm PRIVATE ${CMAKE_SOURCE_DIR}/src/NumericalAlgebra ${CMAKE_SOURCE_DIR}/src/General)
if(arpackng_FOUND)

    target_include_directories(Algorithm PRIVATE ARPACK::ARPACK)
    target_link_libraries(Algorithm PRIVATE ARPACK::ARPACK)

endif()