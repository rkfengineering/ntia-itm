file(GLOB "NTIA_ITM_SOURCES" src/*.cpp)
file(GLOB "NTIA_ITM_HEADERS" include/ITM/*.h)

add_library(ITMLib STATIC ${NTIA_ITM_SOURCES} ${NTIA_ITM_HEADERS})

target_include_directories(ITMLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (NTIA_ITM_BUILD_TESTS)
    add_subdirectory(tests)
endif()