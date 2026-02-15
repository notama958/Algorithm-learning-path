cmake_minimum_required(VERSION 3.24)
project(beacon_program)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  

add_compile_options(-std=c++20)

set(SOURCES
    course_code/mainwindow.cc
    course_code/mainprogram.cc    
    datastructures.cc
    )
    
set(HEADERS
    course_code/mainwindow.hh
    course_code/mainprogram.hh    
    datastructures.hh
    )

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

target_include_directories(${PROJECT_NAME} PRIVATE 
    .
    course_code
)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -pedantic -g -O0)

if(GUI)
# Link to Qt
    find_package(Qt6 QUIET COMPONENTS Core Widgets)
    if(Qt6_FOUND)
        target_link_libraries(${PROJECT_NAME} Qt6::Core Qt6::Widgets)
    else()
        message(STATUS "Qt6 not found")
    endif()

    set_target_properties(${PROJECT_NAME} PROPERTIES
        AUTOMOC ON
        AUTOUIC ON
        AUTORCC ON
    )
endif()