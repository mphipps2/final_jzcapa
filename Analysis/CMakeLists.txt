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

add_executable(AnalysisExample ${ANALYSIS_SRC_DIR}/Containers.cpp ${ANALYSIS_SRC_DIR}/XMLSettingsReader.cpp ${ANALYSIS_SRC_DIR}/DataReader.cpp ${ANALYSIS_SRC_DIR}/WFAnalysis.cpp userFunctions/AnalysisExample.cpp)
TARGET_LINK_LIBRARIES(AnalysisExample ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES(AnalysisExample ${XERCESC_LIBRARY})


install(TARGETS AnalysisExample RUNTIME DESTINATION bin)

