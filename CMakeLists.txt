cmake_minimum_required(VERSION 3.10)
project(velodep VERSION 0.2.0 LANGUAGES CXX)

# Set the CMP0095 policy
cmake_policy(SET CMP0095 NEW)

if(DEFINED ENV{CC})
    set(CMAKE_C_COMPILER $ENV{CC})
endif()
if(DEFINED ENV{CXX})
    set(CMAKE_CXX_COMPILER $ENV{CXX})
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
add_library(velodep SHARED ${SOURCE_FILES})

# If on a Windows platform, define the macro BUILDING_DLL to export APIs
if(WIN32)
    target_compile_definitions(velodep PRIVATE BUILDING_DLL)
endif()


# Set the output directory of binary files
set_target_properties(velodep PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/Debug/bin
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/Release/bin
    LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/Debug/lib
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/Release/lib
    ARCHIVE_OUTPUT_DIRECTORY_DEBUG ${PROJECT_BINARY_DIR}/Debug/lib
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${PROJECT_BINARY_DIR}/Release/lib
)
if(UNIX)
    # Set RPATH (the search path of runtime library)
    set_target_properties(velodep PROPERTIES
        INSTALL_RPATH "\${ORIGIN}"
    )
endif()


# Display the information of compilers
message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "CXX Compiler: ${CMAKE_CXX_COMPILER}")
# Display the information of the generator
message(STATUS "Generator: ${CMAKE_GENERATOR}")
# Display the build type
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")


# Install targets
install(TARGETS velodep
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


# Declare the PACK_DEPENDENCIES option
option(PACK_DEPENDENCIES "Pack required system libraries." OFF)
# Conditionally install the required system libraries if PACK_DEPENDENCIES is defined
if(PACK_DEPENDENCIES)
    # Include the InstallRequiredSystemLibraries module
    include(InstallRequiredSystemLibraries)
    # Debug output to check the contents of CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS
    message(STATUS "CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS: ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}")
    # Install the required system libraries
    if(WIN32)
        install(PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
                DESTINATION bin)
    elseif(UNIX)
        install(PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}
                DESTINATION lib)
    endif()
endif()


# CPack configurations
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
    set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
    set(CPACK_PACKAGE_CONTACT "zhangsupeng18@mails.ucas.ac.cn")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY 
        "velodep - A tool to solve the velocity dependency of the slip behavior 
         and aperture of a fault.")

    if(WIN32)
        set(CPACK_GENERATOR "ZIP")
    elseif(UNIX)
        set(CPACK_GENERATOR "TGZ")
    endif()

    # Include the CPack module
    include(CPack)
 endif()