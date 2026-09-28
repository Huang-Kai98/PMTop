file(GLOB SOURCES "*.cpp")

add_executable(PMTop ${SOURCES})



if (MKL_FOUND)
target_compile_options(PMTop PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(PMTop PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(PMTop PUBLIC $<LINK_ONLY:MKL::MKL>)
endif()


target_include_directories(PMTop PRIVATE  ${CMAKE_SOURCE_DIR}/src/NumericalAlgebra ${CMAKE_SOURCE_DIR}/src/FEM ${CMAKE_SOURCE_DIR}/src/General
        ${CMAKE_SOURCE_DIR}/src/${CMAKE_SOURCE_DIR}/src/Algorithm)
target_link_libraries(PMTop PRIVATE NumericalAlgebra FEM General Algorithm)


if(arpackng_FOUND)

target_include_directories(PMTop PRIVATE ARPACK::ARPACK)
target_link_libraries(PMTop PRIVATE ARPACK::ARPACK) 

endif()

if(nlohmann_json_FOUND)
    target_link_libraries(PMTop PRIVATE nlohmann_json::nlohmann_json)
endif ()

if(spdlog_FOUND)
    target_link_libraries(PMTop PRIVATE spdlog::spdlog)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(PMTop PUBLIC OpenMP::OpenMP_CXX)
endif()