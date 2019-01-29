#Basic Cmake structure created by R.Longo on 12/13/2018
project(JCaPA)
cmake_minimum_required(VERSION 2.8)
#aux_source_directory(. SRC_LIST)

#Colors definition for messages
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

message("${BoldBlue}
      #  #####           ######     #
      # #     #          #     #   # #
      # #                #     #  #   #
      # #        ####    ######  #     #
#     # #       #    #   #       #######
#     # #     # #    #   #       #     #
 #####   #####   ######  #       #     #
${ColourReset}")

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

#To be commented out on old linux...
add_compile_options(-std=c++11)
#add_compile_options(-arch x86_64)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmakeModules/)

#Root support
find_package (ROOT REQUIRED)
include_directories (${ROOT_INCLUDE_DIR})

#Xerces-C support
find_package (Xerces REQUIRED)
include_directories (${XERCESC_INCLUDE})



#Test executable - does nothing
add_executable(test_JCaPA  main.cpp)
TARGET_LINK_LIBRARIES(test_JCaPA ${ROOT_LIBRARIES})
TARGET_LINK_LIBRARIES(test_JCaPA ${XERCESC_LIBRARY})



#Subdirectory structure
#Analysis libraries and user functions 
add_subdirectory(Analysis)

#MonteCarlo folder. The installation will be conditional and disabled by default (not everyone needs to work with it and requires dependencies)
#The starting point will be the work done by Mike Phipps

if (${JCaPA_ENABLE_MC} MATCHES "YES")
    message("${BoldGreen} Monte Carlo extension enabled! ${ColourReset}")
    add_subdirectory(MonteCarlo)

endif (${JCaPA_ENABLE_MC} MATCHES "YES")

#Install test
install(TARGETS test_JCaPA RUNTIME DESTINATION bin)
#Install doxygen files
install(FILES doxygen/JCaPA_doxy.cnf DESTINATION doxygen)

