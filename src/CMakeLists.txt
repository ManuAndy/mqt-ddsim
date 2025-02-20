add_subdirectory("${PROJECT_SOURCE_DIR}/extern/mqt-core" "extern/mqt-core"
                 EXCLUDE_FROM_ALL)

add_library(
  ${PROJECT_NAME}
  ${PROJECT_SOURCE_DIR}/include/Simulator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Simulator.cpp
  ${PROJECT_SOURCE_DIR}/include/CircuitSimulator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CircuitSimulator.cpp
  ${PROJECT_SOURCE_DIR}/include/GroverSimulator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/GroverSimulator.cpp
  ${PROJECT_SOURCE_DIR}/include/ShorSimulator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ShorSimulator.cpp
  ${PROJECT_SOURCE_DIR}/include/ShorFastSimulator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ShorFastSimulator.cpp
  ${PROJECT_SOURCE_DIR}/include/StochasticNoiseSimulator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/StochasticNoiseSimulator.cpp
  ${PROJECT_SOURCE_DIR}/include/DeterministicNoiseSimulator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DeterministicNoiseSimulator.cpp
  ${PROJECT_SOURCE_DIR}/include/HybridSchrodingerFeynmanSimulator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/HybridSchrodingerFeynmanSimulator.cpp
  ${PROJECT_SOURCE_DIR}/include/UnitarySimulator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/UnitarySimulator.cpp
  ${PROJECT_SOURCE_DIR}/include/PathSimulator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/PathSimulator.cpp
  ${PROJECT_SOURCE_DIR}/include/EntanglementSplitSimulator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/EntanglementSplitSimulator.cpp)
target_include_directories(
  ${PROJECT_NAME}
  PUBLIC $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>)

# cmake-lint: disable=C0103
set(JSON_MultipleHeaders
    OFF
    CACHE INTERNAL "")
target_link_libraries(${PROJECT_NAME} PUBLIC MQT::CoreDD)

set(TF_BUILD_TESTS
    OFF
    CACHE BOOL "")
set(TF_BUILD_EXAMPLES
    OFF
    CACHE BOOL "")
set(TF_BUILD_PROFILER
    OFF
    CACHE BOOL "")
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/taskflow" "extern/taskflow")
target_link_libraries(${PROJECT_NAME} PUBLIC Taskflow)
# the following sets the SYSTEM flag for the include dirs of the taskflow libs
# cmake-lint: disable=C0307
set_target_properties(
  Taskflow PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                      $<TARGET_PROPERTY:Taskflow,INTERFACE_INCLUDE_DIRECTORIES>)

add_library(MQT::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
