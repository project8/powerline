cmake_minimum_required (VERSION 2.6)

# Specify the local directory for CMake modules
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

# Define the project
project (powerline)

# Setup the default install prefix
# This gets set to the binary directory upon first configuring.
# If the user changes the prefix, but leaves the flag OFF, then it will remain as the user specified.
# If the user wants to reset the prefix to the default (i.e. the binary directory), then the flag should be set ON.
if (NOT DEFINED SET_INSTALL_PREFIX_TO_DEFAULT)
    set (SET_INSTALL_PREFIX_TO_DEFAULT ON)
endif (NOT DEFINED SET_INSTALL_PREFIX_TO_DEFAULT)
if (SET_INSTALL_PREFIX_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR} CACHE PATH "Install prefix" FORCE)
    set (SET_INSTALL_PREFIX_TO_DEFAULT OFF CACHE BOOL "Reset default install path when when configuring" FORCE)
endif (SET_INSTALL_PREFIX_TO_DEFAULT)

set (LIBS)

# Monarch
set (MONARCH_LOCATION ${CMAKE_SOURCE_DIR}/../katydid/cbuild-develop/Monarch)
find_library (MONARCH_LIBRARIES NAMES MonarchCore PATHS ${MONARCH_LOCATION}/lib DOC "Monarch library")
find_path(MONARCH_INCLUDE_DIR NAMES Monarch.hpp PATHS ${MONARCH_LOCATION}/include DOC "Monarch include directory")
list (APPEND LIBS ${MONARCH_LIBRARIES})
include_directories (${MONARCH_INCLUDE_DIR})

# FFTWf
find_package(FFTWf REQUIRED)
include_directories (${FFTWf_INCLUDE_DIR})
list (APPEND LIBS ${FFTWf_LIBRARIES})
if (FFTWf_THREADS_FOUND)
    set (FFTWf_NTHREADS 1 CACHE STRING "Number of threads to use for FFTWf processes")
    add_definitions (-DFFTWf_NTHREADS=${FFTW_NTHREADS})
    message (STATUS "FFTWf configured to use up to ${FFTWf_NTHREADS} threads.")
else (FFTWf_THREADS_FOUND)
    remove_definitions (-DFFTWf_NTHREADS=${FFTWf_NTHREADS})
endif (FFTWf_THREADS_FOUND)

# Powerline
include_directories (${PROJECT_SOURCE_DIR})

####################################################
# Executables
####################################################

add_executable (powerline powerline.cc)
target_link_libraries (powerline ${LIBS})

add_executable (sweepline sweepline.cc)
target_link_libraries (sweepline ${LIBS})

#add_executable (viewcandidate_elektronjager viewcandidate_elektronjager.cc Waterfall.cc correline_utils.cc Histogram.cc)
#target_link_libraries (viewcandidate_elektronjager ${LIBS})

add_executable (view_candidate view_candidate.cc Waterfall.cc correline_utils.cc)
target_link_libraries (view_candidate ${LIBS})

add_executable (correline correline.cc Waterfall.cc)
target_link_libraries (correline ${LIBS})

add_executable (powerline_twochannel powerline_twochannel.cc Waterfall.cc)
target_link_libraries (powerline_twochannel ${LIBS})

add_executable (dpph_search dpph_search.cc Waterfall.cc)
target_link_libraries (dpph_search ${LIBS})

add_executable (data_simulator data_simulator.cc correline_utils.cc)
target_link_libraries (data_simulator ${LIBS})

add_executable (powerline_elektronjager powerline_elektronjager.cc Waterfall.cc Histogram.cc correline_utils.cc)
target_link_libraries (powerline_elektronjager ${LIBS})

add_executable (correline_electron_huntress correline_electron_huntress.cc Waterfall.cc Histogram.cc correline_utils.cc)
target_link_libraries (correline_electron_huntress ${LIBS})



