# really no optimization in debug mode
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  set(DEBUG_FLAGS "-O0 -Wall -Wextra")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${DEBUG_FLAGS}")
  # suppress warnings in release mode (just not to confuse people)
  set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -w")
else()
  message(FATAL_ERROR "Only GNU Fortran compiler supported")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()

# figure out compile flags
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
set(DEFAULT_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_FLAGS_${BUILD_TYPE}})
