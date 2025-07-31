add_library(crelib_lib STATIC
    crelib.cpp
)

target_include_directories(crelib_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)