set(FasterCap_Sources
    ${CMAKE_CURRENT_LIST_DIR}/AutoConductor.cpp
    ${CMAKE_CURRENT_LIST_DIR}/AutoConductor.h
    ${CMAKE_CURRENT_LIST_DIR}/AutoElement.cpp
    ${CMAKE_CURRENT_LIST_DIR}/AutoElement.h
    ${CMAKE_CURRENT_LIST_DIR}/AutoPanel.cpp
    ${CMAKE_CURRENT_LIST_DIR}/AutoPanel.h
    ${CMAKE_CURRENT_LIST_DIR}/Autorefine.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Autorefine.h
    ${CMAKE_CURRENT_LIST_DIR}/AutoSegment.cpp
    ${CMAKE_CURRENT_LIST_DIR}/AutoSegment.h
    ${CMAKE_CURRENT_LIST_DIR}/MultiplyHierarchical.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MultiplyHierarchical.h
    ${CMAKE_CURRENT_LIST_DIR}/Potential.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Potential.h
    ${CMAKE_CURRENT_LIST_DIR}/SolveCapacitance.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SolveCapacitance.h
    ${CMAKE_CURRENT_LIST_DIR}/SolverGlobal.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SolverGlobal.h
    ${CMAKE_CURRENT_LIST_DIR}/SphereGen.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SphereGen.h
)
   
    
#set(LinAlgebra_Required_Libs
#   debug MSVCRTD.LIB
#   debug MSVCPRTD.LIB
#   optimized MSVCRT.LIB
#   optimized MSVCPRT.LIB)

# not creating a library (to save on compile time)
#add_library(Solver ${Solver_Sources})
# libraries that needs to be linked to the present library
#target_link_libraries(Solver ${Solver_Required_Libs})
#target_link_libraries(Solver ${EXTRA_LIBS})


