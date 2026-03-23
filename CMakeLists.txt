cmake_minimum_required(VERSION 3.16)
project(GAGG_Calorimeter)


find_package(Geant4 REQUIRED ui_all vis_all)


include(${Geant4_USE_FILE})


include_directories(${PROJECT_SOURCE_DIR}/include)


file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


file(COPY ${PROJECT_SOURCE_DIR}/macros DESTINATION ${PROJECT_BINARY_DIR})


add_executable(GAGG_sim main.cc ${sources} ${headers})
target_link_libraries(GAGG_sim ${Geant4_LIBRARIES})
