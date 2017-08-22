PROJECT(stap++)

cmake_minimum_required(VERSION 3.0)

OPTION(STAP++_DEBUG "Output detailed results for debugging." OFF)
IF(STAP++_DEBUG)
   ADD_DEFINITIONS(-D_DEBUG_)
ENDIF()

INCLUDE_DIRECTORIES(h)

AUX_SOURCE_DIRECTORY(cpp SRC)
FILE(GLOB_RECURSE HEAD h/*.h)

source_group(SOURCE\ FILES FILES ${SRC})
source_group(HEADER\ FILES FILES ${HEAD})

ADD_EXECUTABLE(stap++ ${SRC} ${HEAD})
