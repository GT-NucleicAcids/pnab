# Try to find OpenBabel3 headers and libraries
# Defines:
#  OPENBABEL3_FOUND - system has OpenBabel3
#  OPENBABEL3_INCLUDE_DIR - the OpenBabel3 include directory
#  OPENBABEL3_LIBRARY - Link these to use OpenBabel3
#  IF OPENBABEL_DIR is defined, will look there first

if(OPENBABEL3_INCLUDE_DIR AND OPENBABEL3_LIBRARY)
  # in cache already or user-specified
  set(OPENBABEL3_FOUND TRUE)

else()

  if(NOT OPENBABEL3_INCLUDE_DIR)
      find_path(OPENBABEL3_INCLUDE_DIR openbabel/obconversion.h
        PATHS
        ${OPENBABEL_DIR}/include/openbabel3
        ${OPENBABEL_DIR}/include
        ${OPENBABEL_DIR}/Library/include/openbabel3
      )
    if(OPENBABEL3_INCLUDE_DIR)
      message(STATUS "Found Open Babel include files at ${OPENBABEL3_INCLUDE_DIR}")
      set(OPENBABEL3_INCLUDE_DIR ${OPENBABEL3_INCLUDE_DIR})
    endif()
  endif()

  if(NOT OPENBABEL3_LIBRARY)
  find_library(OPENBABEL3_LIBRARY NAMES openbabel openbabel3 openbabel-3
      PATHS
      ${OPENBABEL_DIR}/lib
      ${OPENBABEL_DIR}/Library/bin
    )
    if(OPENBABEL3_LIBRARY)
      message(STATUS "Found Open Babel library at ${OPENBABEL3_LIBRARY}")
      set(OPENBABEL3_LIBRARY ${OPENBABEL3_LIBRARY})
    endif()
  endif()

  if(OPENBABEL3_INCLUDE_DIR AND OPENBABEL3_LIBRARY)
    set(OPENBABEL3_FOUND TRUE)
  endif()

  mark_as_advanced(OPENBABEL3_INCLUDE_DIR OPENBABEL3_LIBRARY)
endif()
