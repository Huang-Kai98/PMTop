file(GLOB SOURCES "*.cpp")



add_library(NumericalAlgebra ${SOURCES}
        ComplexScalar.hpp
        ComplexVector.hpp
        ComplexMatrix.cpp)
target_include_directories(NumericalAlgebra PRIVATE ${CMAKE_SOURCE_DIR}/src/General)
target_link_libraries(NumericalAlgebra PRIVATE General) 

if (MKL_FOUND)
include_directories(${MKL_INCLUDE_DIRS})
target_compile_options(NumericalAlgebra PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(NumericalAlgebra PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(NumericalAlgebra PUBLIC $<LINK_ONLY:MKL::MKL>)
endif()

if(arpackng_FOUND)
target_include_directories(NumericalAlgebra PRIVATE ARPACK::ARPACK)
target_link_libraries(NumericalAlgebra PRIVATE ARPACK::ARPACK)
endif()
