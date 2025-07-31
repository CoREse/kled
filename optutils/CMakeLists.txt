add_library(optutils_lib STATIC
    OptHelper.cpp
)

target_include_directories(optutils_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)