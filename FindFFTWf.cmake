# Find the FFTWf includes and library.
# 
# This module defines
# FFTWf_INCLUDE_DIR, where to locate fftw3.h file
# FFTWf_LIBRARIES, the libraries to link against to use fftw3
# FFTWf_FOUND.  If false, you cannot build anything that requires fftw3.
# FFTWf_LIBRARY, where to find the libfftw3 library.
# FFTWf_THREADS_FOUND, If true, can use multithreaded FFTWf routines.
# FFTWf_THREADS_LIBRARY, where to find the libfftw3_threads library.

set(FFTWf_FOUND FALSE)
set(FFTWf_THREADS_FOUND FALSE)
if(FFTWf_LIBRARY AND FFTWf_INCLUDE_DIR)
  set(FFTWf_FIND_QUIETLY TRUE)
endif()

find_path(FFTWf_INCLUDE_DIR fftw3.h
  $ENV{FFTW_DIR}/include
  $ENV{FFTW3} $ENV{FFTW3}/include $ENV{FFTW3}/api
  /usr/local/include
  /usr/include
  /opt/fftw3/include
  DOC "Specify the directory containing fftw3.h"
)

find_library(FFTWf_LIBRARY NAMES fftw3f fftw3-3f PATHS
  $ENV{FFTW_DIR}/lib
  $ENV{FFTW3} $ENV{FFTW3}/lib $ENV{FFTW3}/.libs
  /usr/local/lib
  /usr/lib 
  /opt/local/lib
  /opt/fftw3/lib
  DOC "Specify the fftw3f library here."
)

find_library(FFTWf_THREADS_LIBRARY NAMES fftw3f_threads fftw3-3f_threads PATHS
  $ENV{FFTW_DIR}/lib
  $ENV{FFTW3} $ENV{FFTW3}/lib $ENV{FFTW3}/.libs
  /usr/local/lib
  /usr/lib 
  /opt/local/lib
  /opt/fftw3/lib
  DOC "Specify the fftw3f threads library here."
)

# make sure pthreads is present, and add to list of threads libraries
if(FFTWf_THREADS_LIBRARY)
  find_package(Threads)
  if(CMAKE_USE_PTHREADS_INIT)
    set(FFTWf_THREADS_LIBRARY ${FFTWf_THREADS_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
  else()
    message(WARNING "Found FFTWf threads but no pthreads! Will not use FFTWf threads.")
    unset(FFTWf_THREADS_LIBRARY)
  endif()
endif()


if(FFTWf_INCLUDE_DIR AND FFTWf_LIBRARY)
  set(FFTWf_FOUND TRUE)
  if(NOT FFTWf_FIND_QUIETLY)
     message(STATUS "Found fftw3f includes at ${FFTWf_INCLUDE_DIR}")
     message(STATUS "Found fftw3f library at ${FFTWf_LIBRARY}")
  endif()
  
  if(FFTWf_THREADS_LIBRARY)
    set(FFTWf_THREADS_FOUND TRUE)
    if(NOT FFTWf_FIND_QUIETLY)
      message(STATUS "Found fftw3f threads and pthreads libraries at ${FFTWf_THREADS_LIBRARY}")
    endif()
  endif()
endif()

mark_as_advanced(FFTWf_FOUND FFTWf_LIBRARY FFTWf_INCLUDE_DIR FFTWf_THREADS_FOUND)
if (FFTWf_THREADS_LIBRARY)
  mark_as_advanced(FFTWf_THREADS_LIBRARY)
  set(FFTWf_LIBRARIES ${FFTWf_LIBRARY} ${FFTWf_THREADS_LIBRARY})
else()
  set(FFTWf_LIBRARIES ${FFTWf_LIBRARY})
endif()
