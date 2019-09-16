if(${CMAKE_SYSTEM_NAME} MATCHES "Windows") # we'll assume 'Windows' means MSVC for simplicity
        # boost links via preprocessor directives, so we need to set the dir in which the libs are located.
        link_directories(${USD_LIB_DIR})
endif()

add_executable (squishCube
        main.cpp
        )

find_package(PythonLibs)

target_include_directories(squishCube
        PUBLIC
            ${USD_INC_DIR}
            ${PYTHON_INCLUDE_PATH}
            ${USD_BOOST_INC_DIR}
        )

 set (LINK_LIBS 
                ${USD_LIB_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gf${LIB_SUFFIX}
                ${USD_LIB_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}sdf${LIB_SUFFIX}
                ${USD_LIB_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}tf${LIB_SUFFIX}
                ${USD_LIB_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}usd${LIB_SUFFIX}
                ${USD_LIB_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}usdGeom${LIB_SUFFIX}
                ${USD_LIB_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}vt${LIB_SUFFIX}
                ${PYTHON_LIBRARY}
)

if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows") 
        # this library is not needed on windows builds (linked via pragmas aready by boost)
        list(APPEND LINK_LIBS ${USD_LIB_DIR}/libboost_python.so) # todo find library (but appears to not be necessary)
endif()
       

target_link_libraries(squishCube ${LINK_LIBS})
