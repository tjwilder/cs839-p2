cmake_minimum_required (VERSION 2.8.11)
set (CMAKE_CXX_STANDARD 11)
set(Python_ADDITIONAL_VERSIONS 2.7)
find_package(PythonInterp)

project (USD-Tests)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows") # we'll assume 'Windows' means MSVC for simplicity
    # windows.h defines max as a macro, which interferes with std::limits<>::max()
    add_compile_definitions(NOMINMAX)
    # disable warnings about some casting inconsistencies
    set (CMAKE_CXX_FLAGS "/wd\"4244\" /wd\"4305\" /wd\"4267\" ${CMAKE_CXX_FLAGS}")
    # lib suffix is different on windows than ubuntu
    # surely cmake can take care of this automatically, but I don't know how.
    set (LIB_SUFFIX ".lib")
else()
    set (CMAKE_CXX_FLAGS "-Wno-deprecated ${CMAKE_CXX_FLAGS}")
    set (LIB_SUFFIX ".so")
endif()

add_definitions(-DBUILD_COMPONENT_SRC_PREFIX="" -DBUILD_OPTLEVEL_DEV)

set(USD_ROOT_DIR "K:/Program Files/USD")
set(USD_LIB_DIR ${USD_ROOT_DIR}/lib/)
set(USD_INC_DIR ${USD_ROOT_DIR}/include/)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows") # we'll assume 'Windows' means MSVC for simplicity
    set(USD_BOOST_INC_DIR ${USD_INC_DIR}/boost-1_65_1)
endif()

add_subdirectory(tests)
