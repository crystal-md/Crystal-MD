# this file is included in src directory.

set(CURRENT_ARCH_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/arch)

# set ARCH_NAME and config header file
include(${CURRENT_ARCH_SOURCE_DIR}/arch_configure.cmake)

# set ARCH_SRC_PATH
if (MD_SUNWAY_ARCH_ENABLE_FLAG) #  Sunway
    if (NOT MD_SUNWAY_ARCH_SRC_PATH STREQUAL "")
        set(ARCH_SRC_PATH ${MD_SUNWAY_ARCH_SRC_PATH})
    endif ()
elseif (MD_TH_ARCH_ENABLE_FLAG) # TH-3
    if (NOT MD_TH_ARCH_SRC_PATH STREQUAL "")
        set(ARCH_SRC_PATH ${MD_TH_ARCH_SRC_PATH})
    endif ()
elseif (MD_CUDA_ARCH_ENABLE_FLAG) # CUDA
    if (NOT MD_CUDA_ARCH_SRC_PATH STREQUAL "")
        set(ARCH_SRC_PATH ${MD_CUDA_ARCH_SRC_PATH})
    endif ()
elseif (MD_HIP_ARCH_ENABLE_FLAG) # HIP
    if (NOT MD_HIP_ARCH_SRC_PATH STREQUAL "")
        set(ARCH_SRC_PATH ${MD_HIP_ARCH_SRC_PATH})
    endif ()
endif ()

# set ARCH_LIBS and default ARCH_SRC_PATH
if (ARCH_NAME)
    set(ARCH_LIBS md_arch_${ARCH_NAME})  # may link the arch lib to libmd
    # set default arch directory
    if (NOT ARCH_SRC_PATH)
        set(ARCH_SRC_PATH "${CURRENT_ARCH_SOURCE_DIR}/${ARCH_NAME}")
    endif ()
    # check ARCH_SRC_PATH and add_subdirectory
    if (NOT EXISTS "${ARCH_SRC_PATH}" OR NOT IS_DIRECTORY "${ARCH_SRC_PATH}")
        message(FATAL_ERROR "Architecture source files directory not found: ${ARCH_SRC_PATH}")
    else ()
        MESSAGE(STATUS "Arch source files is ${ARCH_SRC_PATH}")
#        set(MD_SOURCE_INCLUDES "${PROJECT_SOURCE_DIR}/src" CACHE PATH "PATH of includes in arch code.")
        add_subdirectory(${ARCH_SRC_PATH} arch_${ARCH_NAME})
    endif ()
endif ()

# ARCH FILES
set(ARCH_FILES
        ${CURRENT_ARCH_SOURCE_DIR}/arch_building_config.h
        ${CURRENT_ARCH_SOURCE_DIR}/arch_env.hpp
        ${CURRENT_ARCH_SOURCE_DIR}/arch_imp.h
        ${CURRENT_ARCH_SOURCE_DIR}/arch_macros.h
        ${CURRENT_ARCH_SOURCE_DIR}/hardware_accelerate.hpp
        )
