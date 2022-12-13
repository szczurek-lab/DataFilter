cmake_minimum_required(VERSION 3.10)

project(datafilter VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)

# Whether to use conan to install dependencies: ON for yes, OFF for no
set(USE_CONAN ON CACHE BOOL "Whether to use conan to install dependencies: ON (default) or OFF")

# For debugging: memory leakage
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")

# Set directories.
set(SRC_DIR ${CMAKE_SOURCE_DIR}/source)
set(HEADER_DIR ${CMAKE_SOURCE_DIR}/header)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Required libraries.
if(USE_CONAN)
    include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
    conan_basic_setup()
else()
    find_package(Boost COMPONENTS program_options REQUIRED)
    find_package(dlib REQUIRED)
    find_package(ZLIB REQUIRED)
    find_package(Threads REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS})
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

# Include header files of this project.
include_directories(${HEADER_DIR})

# Generate library.
add_executable(
        ${PROJECT_NAME}
        ${SRC_DIR}/main.cpp
        ${SRC_DIR}/process_data.cpp
        ${SRC_DIR}/config.cpp
        ${SRC_DIR}/output_data.cpp
        ${SRC_DIR}/read_counts_file.cpp
        ${SRC_DIR}/probabilities.cpp
        )
target_compile_options(${PROJECT_NAME} PRIVATE -fPIC)

if(USE_CONAN)
    target_link_libraries(${PROJECT_NAME} PRIVATE ${CONAN_LIBS})
else()
    target_link_libraries(${PROJECT_NAME} PRIVATE dlib::dlib ZLIB::ZLIB ${Boost_LIBRARIES})
endif()

if (UNIX AND NOT APPLE)
    target_link_libraries(${PROJECT_NAME} PRIVATE -lm -ldl)
endif ()
