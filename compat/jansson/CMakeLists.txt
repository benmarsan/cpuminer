# jansson library
set(jansson_SOURCES
        config.h
        dump.c
        hashtable.c
        hashtable.h
        jansson.h
        jansson_private.h
        load.c
        strbuffer.c
        strbuffer.h
        utf.c
        utf.h
        util.h
        value.c
)

add_library(jansson ${jansson_SOURCES})
target_include_directories(jansson PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
