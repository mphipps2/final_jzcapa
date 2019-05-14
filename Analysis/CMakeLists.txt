message("${BoldYellow}----------------------------")
message("Preparing Analysis libraries...")
message("----------------------------${ColourReset}")

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Analysis/include)

#In the following the list of source files to be linked
set (ANALYSIS_SRC_DIR "${CMAKE_SOURCE_DIR}/Analysis/src")

#### Get all header files
FILE(GLOB_RECURSE PA_INC
    "${CMAKE_SOURCE_DIR}/Analysis/include/*.hh"
)
FILE(GLOB_RECURSE DEPENDENCIES
    "${CMAKE_SOURCE_DIR}/Analysis/include/*.hh"
    "${CMAKE_SOURCE_DIR}/Analysis/src/*.cpp"
)

add_executable(AnalysisExample ${DEPENDENCIES} userFunctions/AnalysisExample.cpp)
TARGET_LINK_LIBRARIES(AnalysisExample ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES(AnalysisExample ${XERCESC_LIBRARY})

install(TARGETS AnalysisExample RUNTIME DESTINATION bin)


add_executable(myAnalysis ${DEPENDENCIES} userFunctions/myAnalysis.cpp)
TARGET_LINK_LIBRARIES(myAnalysis ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES(myAnalysis ${XERCESC_LIBRARY})

install(TARGETS myAnalysis RUNTIME DESTINATION bin)


add_executable(postAnalysis ${DEPENDENCIES} userFunctions/postAnalysis.cpp)
TARGET_LINK_LIBRARIES(postAnalysis ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES(postAnalysis ${XERCESC_LIBRARY})

install(TARGETS postAnalysis RUNTIME DESTINATION bin)
