cmake_minimum_required(VERSION 3.10)

project(datafilter VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)

# Whether to use conan to install dependencies: ON for yes, OFF for no
set(USE_CONAN OFF CACHE BOOL "Whether to use conan to install dependencies: ON or OFF (default)")

# For debugging: memory leakage
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")

# Set directories.
include(GNUInstallDirs)
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(INC_DIR ${CMAKE_SOURCE_DIR}/include/datafilter)
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
include_directories(${INC_DIR})

# Generate library.
add_library(
        ${PROJECT_NAME}
        SHARED
        ${SRC_DIR}/process_data.cpp
        ${SRC_DIR}/config.cpp
        ${SRC_DIR}/output_data.cpp
        ${SRC_DIR}/read_counts_file.cpp
        ${SRC_DIR}/probabilities.cpp
)
target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

if(USE_CONAN)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${CONAN_LIBS})
else()
    target_link_libraries(${PROJECT_NAME} PUBLIC dlib::dlib ZLIB::ZLIB ${Boost_LIBRARIES})
endif()

if (UNIX AND NOT APPLE)
    target_link_libraries(${PROJECT_NAME} PUBLIC -lm -ldl)
endif ()

# Install library.
install(
        TARGETS
        ${PROJECT_NAME}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Install include directory.
install(
        DIRECTORY
        ${INC_DIR}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Generate executable.
add_executable(
        ${PROJECT_NAME}_bin
        ${SRC_DIR}/main.cpp
)
target_link_libraries(${PROJECT_NAME}_bin PRIVATE ${PROJECT_NAME})
set_target_properties(
        ${PROJECT_NAME}_bin
        PROPERTIES
        OUTPUT_NAME ${PROJECT_NAME}
        BUILD_WITH_INSTALL_RPATH TRUE
        INSTALL_RPATH_USE_LINK_PATH TRUE
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}:${INSTALL_RPATH}"
)

# Install executable.
install(
        TARGETS
        ${PROJECT_NAME}_bin
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
