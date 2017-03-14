cmake_minimum_required(VERSION 2.6)
project(ppcg.exe)

macro(use_cxx11)
    if (CMAKE_VERSION VERSION_LESS "3.1")
        if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            set (CMAKE_CXX_FLAGS "--std=gnu++11 -fpermissive -O3 ${CMAKE_CXX_FLAGS}")
        endif ()
    else ()
        set (CMAKE_CXX_STANDARD 11)
    endif ()
endmacro(use_cxx11)


set(CMAKE_CXX_STANDARD 11)
use_cxx11()


set(SOURCE_FILES
    UncollapsedSampler/igmm.cpp
    UncollapsedSampler/Restaurant.cpp
    UncollapsedSampler/Restaurant.h
    UncollapsedSampler/Table.cpp
    UncollapsedSampler/Table.h)

INCLUDE_DIRECTORIES(../GMMBase/)
LINK_DIRECTORIES(../GMMBase)
INCLUDE_DIRECTORIES(../FastMat/)
LINK_DIRECTORIES(../FastMat)
find_package (Threads)
add_executable(ppcg.exe ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(ppcg.exe FastMat)
TARGET_LINK_LIBRARIES (ppcg.exe ${CMAKE_THREAD_LIBS_INIT})
TARGET_LINK_LIBRARIES(ppcg.exe GMMBase)
