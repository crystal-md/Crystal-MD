#############
## options ##
#############
# change to mpicc and mpicxx
#set(CMAKE_C_COMPILER mpicc -cc=gcc -cxx=g++)
#set(CMAKE_CXX_COMPILER mpicxx -cc=gcc -cxx=g++)

option(OpenMP_ENABLE_FLAG "Use OpenMP" OFF) #change this flag to OFF to disable OpenMP
option(MPI_ENABLE_FLAG "Use MPI library" ON) #change this flag to false to disable mpi
option(TEST_ENABLE_FLAG "Enable test" ON) # enable test
option(TEST_MPI_ENABLE_FLAG "Enable MPI in test" ON) # enable mpi in test, its value depends on option MPI_ENABLE_FLAG.
option(TOOLS_BUILD_ENABLE_FLAG "Enable tools building" ON) # enable tools building (in tools directory) binary.(tools example: convert simulation result binary file to text file)

## architecture ralated values (only used in src/arch dir).
option(SUNWAY_ARCH_ENABLE_FLAG "Enable sunway athread" OFF) # enable sunway athread if its running on sunway system.
option(CUDA_ARCH_ENABLE_FLAG "Enable GPU or DCU using CUDA and HIP" OFF) # enable GPU if its running on nvidia GPU.
set(SUNWAY_ARCH_SRC_PATH "" CACHE STRING "Source files directory of sunway architecture") # source file directory of sunway arch code.
set(CUDA_ARCH_SRC_PATH "" CACHE STRING "Source directory of CUDA/HIP architecture code") # source file directory of cuda arch code.

set(MD_RAND "MT" CACHE STRING "random number generating algorithm") # random number generating
# options are:
#   LCG: linear congruential
#   MT: mersenne twister
#   STC: subtract with carry
#   xoshiro: http://xoshiro.di.unimi.it
#   LEGACY: lammps legacy
#   REAL: real random number privided by linux OS.

if (CMAKE_BUILD_TYPE MATCHES "^(Debug|DEBUG|debug)$")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")

    set(MD_DEV_MODE ON)
endif ()


#############
## const ##
#############
# all variables here start with "MD_"
set(EXECUTE_BIN_NAME CrystalMD)
set(MD_LIB_NAME md) # use PARENT_SCOPE to modify globle variable.
set(MD_FRONTEND_LIB_NAME frontend)

# test
set(MD_TEST_NAME "md-test")
