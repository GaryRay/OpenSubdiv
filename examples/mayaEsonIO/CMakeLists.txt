#
#   Copyright

set(MAYA_FIND_QUIETLY TRUE)

if(NOT MAYA_FOUND)
    message(STATUS
        "not maya found"
    )
    return()
endif()

set(PLATFORM_LIBRARIES
    "${OSD_LINK_TARGET}"
)

include_directories(
    "${PROJECT_SOURCE_DIR}/opensubdiv"
    "${MAYA_INCLUDE_DIRS}"
	"../osd2eson"
)

set(SOURCE_FILES
	main.cpp
	mayaEsonReader.cpp
	mayaEsonWriter.cpp
    mayaEsonIO.cpp
    ../raytrace/convert_bezier.cpp
    ../osd2eson/eson.cc

)

set(HEADER_FILES
	mayaEsonReader.h
	mayaEsonWriter.h
    mayaEsonIO.h
    ../osd2eson/eson.h
)

if(UNIX)
    set(PLATFORM_COMPILE_FLAGS
        -D_BOOL
        -DREQUIRE_IOSTREAM
        -DLINUX
    )
    set(PLATFORM_PLUGIN_EXTENSION
        .so
    )
    set(PLATFORM_LINK_FLAGS
    )
elseif(WIN32)
    set(PLATFORM_COMPILE_FLAGS
        /D_AFXDLL
        /DNT_PLUGIN
        /DREQUIRE_IOSTREAM
    )
    set(PLATFORM_PLUGIN_EXTENSION
        .mll
    )
    set(PLATFORM_LINK_FLAGS
        "/export:initializePlugin /export:uninitializePlugin"
    )
endif()

add_definitions(
    ${PLATFORM_COMPILE_FLAGS}
)

add_library(mayaEsonIO SHARED
    ${SOURCE_FILES}
    ${HEADER_FILES}
)

set_target_properties(mayaEsonIO
    PROPERTIES
    OUTPUT_NAME "mayaEsonIO"
    PREFIX ""
    SUFFIX ${PLATFORM_PLUGIN_EXTENSION}
    LINK_FLAGS "${PLATFORM_LINK_FLAGS}"
)

target_link_libraries(mayaEsonIO
    "${MAYA_Foundation_LIBRARY}"
    "${MAYA_OpenMaya_LIBRARY}"
    "${MAYA_OpenMayaRender_LIBRARY}"
    "${MAYA_OpenMayaUI_LIBRARY}"
    ${PLATFORM_LIBRARIES}
)
