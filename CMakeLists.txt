cmake_minimum_required(VERSION 3.10)
project(rsfsolver VERSION 0.1.1 LANGUAGES CXX)

# Set the compiler based on the system environment
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CMAKE_C_COMPILER /usr/bin/gcc)
    set(CMAKE_CXX_COMPILER /usr/bin/g++)
elseif (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(CMAKE_C_COMPILER cl)
    set(CMAKE_CXX_COMPILER cl)
endif()

set(CMAKE_CXX_STANDARD 20)

# Statically link runtime libraries（optional）
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")

# Set the directory of the included header files
include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB_RECURSE HEADER_H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB_RECURSE HEADER_HPP_FILES ${PROJECT_SOURCE_DIR}/include/*.hpp)
set(HEADER_FILES ${HEADER_H_FILES} ${HEADER_HPP_FILES})
# Set the directory of the source files
file(GLOB_RECURSE SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Creat the dynamic link library
add_library(rsf SHARED ${SOURCE_FILES} ${HEADER_FILES})

# If on a Windows platform, define the macro BUILDING_DLL to export APIs
if (WIN32)
    target_compile_definitions(rsf PRIVATE BUILDING_DLL)
endif()

# # Set the output directory of binary files
set_target_properties(rsf PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/bin/Debug
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/bin/Release
    LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/lib/Debug
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/lib/Release
    ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/lib/Debug
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/lib/Release
)

# Display the information of compilers
message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "CXX Compiler: ${CMAKE_CXX_COMPILER}")
# Display the information of the generator
message(STATUS "Generator: ${CMAKE_GENERATOR}")
# Display the build type
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

# Install targets
install(TARGETS rsf
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# Install header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

# # Install the *.py files in the pylthon directory
# file(GLOB_RECURSE PYTHON_SUB_FILES ${PROJECT_SOURCE_DIR}/python/*.py)
# install(FILES ${PYTHON_SUB_FILES} DESTINATION python)

# Install *.py files to the root directory
file(GLOB_RECURSE PYTHON_FILES ${PROJECT_SOURCE_DIR}/*.py)
# install(FILES ${PROJECT_SOURCE_DIR}/*.py DESTINATION .)
install(FILES ${PYTHON_FILES} DESTINATION .)

# CPack configurations
set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_CONTACT "zhangsupeng18@mails.ucas.ac.cn")
set(CPACK_GENERATOR "ZIP")

# Include the CPack module
include(CPack)