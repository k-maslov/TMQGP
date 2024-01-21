cmake_minimum_required(VERSION 3.1)
project(TMQGP)

include(GenerateExportHeader)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
set(CMAKE_CXX_STANDARD 14)

if (UNIX)
set(SWIG_EXEC swig)
endif(UNIX)

find_package(Python3 COMPONENTS Interpreter Development NumPy)

find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

#IF(MPI_CXX_FOUND)
#        INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
#        LIST(APPEND SCR_EXTERNAL_LIBS ${MPI_CXX_LIBRARIES})
#ENDIF(MPI_CXX_FOUND)

if(WIN32)
#set(SWIG_EXEC "C:\\msys64\\usr\\bin\\swig.exe"
set(SWIG_EXEC "C:\\Users\\const.maslov\\Numerics\\swigwin-4.0.1\\swig.exe")
endif(WIN32) 

if(UNIX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -g -Ofast -std=c++2a -fopenmp")
endif(UNIX)

if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zi /MP3")
endif()

set(VERBOSE TRUE)

file (GLOB_RECURSE libTMQGP_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/TMQGP/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TMQGP/*.h
)

file (GLOB_RECURSE libTMQGP_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/TMQGP/*.h
)

if (UNIX)
# include_directories("/usr/include/levmar")
include_directories(${Python3_INCLUDE_DIRS})
include_directories(${Python3_NumPy_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/TMQGP)
# include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../eosstudio)

link_directories("/usr/lib/levmar")
link_directories("/usr/lib/x86_64-linux-gnu/")
# link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../eosstudio/lib)

endif (UNIX)

if (WIN32)
    include_directories("C:\\Users\\const.maslov\\Numerics\\levmar\\levmar-2.6")
    include_directories("C:\\Users\\const.maslov\\Numerics\\gsl")
    include_directories("C:\\Users\\const.maslov\\Numerics\\K\\c\\mygsl")
    include_directories("C:\\Users\\const.maslov\\Numerics\\K\\include")
    include_directories("C:\\ProgramData\\Anaconda3\\include")
    include_directories("C:\\ProgramData\\Anaconda3\\lib\\site-packages\\numpy\\core\\include")
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/EosStudio)

    include_directories("C:\\Users\\const.maslov\\Numerics\\EosStudio\\packages\\stacktrace.0.1\\build\\native\\include")
    link_directories("C:\\Users\\const.maslov\\Numerics\\gsl\\lib\\x64\\Release")
    link_directories("C:\\Users\\const.maslov\\Numerics\\levmar\\levmar-2.6\\x64\\Release")
    link_directories("C:\\Users\\const.maslov\\Numerics\\lapack")
    link_directories("C:\\ProgramData\\Anaconda3\\libs")
    link_directories("C:\\Users\\const.maslov\\Numerics\\libf2c")
    link_directories("C:\\Users\\const.maslov\\Numerics\\K\\c\\gslbase\\x64\\Debug")
    link_directories("C:\\Users\\const.maslov\\Numerics\\K\\c\\Tdbase\\x64\\Debug")
    link_directories("C:\\Users\\const.maslov\\Numerics\\EosStudio\\packages\\stacktrace.0.1\\build\\native\\lib\\x64\\Release")
endif (WIN32)

set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/lib)

add_library(libTMQGP SHARED ${libTMQGP_SRC})
set_target_properties(libTMQGP PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

#set_property(TARGET EosStudio
#    PROPERTY VS_PACKAGE_REFERENCES "boost_1.7.0"
#)



add_library(TMQGP SHARED ${CMAKE_CURRENT_SOURCE_DIR}/wrap.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/wrap.i)


add_executable(DirectUse ${CMAKE_CURRENT_SOURCE_DIR}/DirectUse.cpp)
# add_executable(tmqgp_mpi ${CMAKE_CURRENT_SOURCE_DIR}/tmqgp_mpi.cpp)


set_target_properties(TMQGP PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

set_target_properties(TMQGP PROPERTIES RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_CURRENT_SOURCE_DIR}/bin")

set_target_properties(TMQGP PROPERTIES PREFIX _)
set_target_properties(TMQGP PROPERTIES IMPORT_PREFIX _)
#set_target_properties(TMQGP PROPERTIES RELW "_rd")
if (WIN32)
set_target_properties(TMQGP PROPERTIES SUFFIX .pyd)
endif (WIN32)
set_target_properties(TMQGP PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
add_dependencies(TMQGP libTMQGP)
add_dependencies(DirectUse libTMQGP)



#set_property(TARGET EosStudio
#    PROPERTY VS_PACKAGE_REFERENCES "stacktrace_0.1"
#)

#set_property(TARGET EosStudio
#    PROPERTY VS_PACKAGE_REFERENCES "stacktrace_0.1"
#)

if (UNIX)
target_link_libraries(TMQGP gsl gslcblas lapack ${Python3_LIBRARIES} libTMQGP ${MPI_CXX_LIBRARIES})
target_link_libraries(libTMQGP gsl gslcblas lapack ${Python3_LIBRARIES} ${MPI_CXX_LIBRARIES})
target_link_libraries(DirectUse gsl gslcblas lapack ${Python3_LIBRARIES} libTMQGP ${MPI_CXX_LIBRARIES})
# target_link_libraries(tmqgp_mpi gsl gslcblas lapack ${Python3_LIBRARIES} libTMQGP ${MPI_CXX_LIBRARIES})
endif(UNIX)


if (WIN32)
target_link_libraries(TMQGP gsl cblas blas blas_nowrap levmar lapack python38 vcf2c libTMQGP gslbase Tdbase StackTrace)
target_link_libraries(libTMQGP gsl cblas blas blas_nowrap levmar lapack python38 vcf2c gslbase Tdbase StackTrace)
endif(WIN32)

# WES stands for wrapEosStudio
get_target_property(WES_RUNTIME_OUTPUT_DIRECTORY TMQGP RUNTIME_OUTPUT_DIRECTORY)

if(WIN32)
link_directories(${LIBRARY_OUTPUT_PATH})
link_directories(${libTMQGP_RUNTIME_OUTPUT_DIRECTORY})
endif(WIN32)

#set (SWIGCOMMAND "${SWIG_EXEC} -DSWIG_PYTHON_LEGACY_BOOL -I${CMAKE_CURRENT_SOURCE_DIR}/EosStudio -c++ -python -o ${CMAKE_CURRENT_SOURCE_DIR}/wrap.cpp 
#-outdir ${EosStudio_RUNTIME_OUTPUT_DIRECTORY}/ ${CMAKE_CURRENT_SOURCE_DIR}/wrap.i" )
set(CMAKE_VERBOSE_MAKEFILE on)
#-I${CMAKE_CURRENT_SOURCE_DIR}/EosStudio 
#-IC:/Users/const.maslov/Numerics/EosStudio/packages/stacktrace.0.1/build/native/include/


add_custom_command(
OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/wrap.cpp
	#message(WARNING "Swig_exec=${SWIG_EXEC}")
COMMAND ${SWIG_EXEC} -DSWIG_PYTHON_LEGACY_BOOL -I${CMAKE_CURRENT_SOURCE_DIR}/EosStudio -c++ -python -o ${CMAKE_CURRENT_SOURCE_DIR}/wrap.cpp -outdir ${CMAKE_CURRENT_SOURCE_DIR}/lib/ ${CMAKE_CURRENT_SOURCE_DIR}/wrap.i
DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/wrap.i
DEPENDS ${libTMQGP_HEADERS}
)

#COMMAND ${SWIG_EXEC} -DSWIG_PYTHON_LEGACY_BOOL -I${GOVNODIR} -I${CMAKE_CURRENT_SOURCE_DIR}/EosStudio -c++ -python -o ${CMAKE_CURRENT_SOURCE_DIR}/wrap.cpp 
 #       -outdir ${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_CFG_INTDIR} ${CMAKE_CURRENT_SOURCE_DIR}/wrap.i

# add_executable(DirectUse DirectUse.cpp)

# if (UNIX)
# target_link_libraries(DirectUse gsl gslcblas levmar lapack EosStudio python3.10)
# endif(UNIX)

# if(WIN32)
# target_link_libraries(DirectUse gsl cblas blas blas_nowrap levmar lapack libTMQGP python38 vcf2c gslbase Tdbase)
# endif(WIN32)
