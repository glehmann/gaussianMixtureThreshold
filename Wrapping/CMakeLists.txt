FIND_PACKAGE(WrapITK REQUIRED)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

BEGIN_WRAPPER_LIBRARY("${PROJECT_NAME}")
SET(WrapITK_VERSION "${WrapITK_VERSION_MAJOR}.${WrapITK_VERSION_MINOR}.${WrapITK_VERSION_PATCH}")
IF("${WrapITK_VERSION}" STRGREATER "0.3")
  SET(WRAPPER_LIBRARY_DEPENDS VXLNumerics Base Numerics)
ELSE("${WrapITK_VERSION}" STRGREATER "0.3")
  SET(WRAPPER_LIBRARY_DEPENDS Base BaseNumerics)
ENDIF("${WrapITK_VERSION}" STRGREATER "0.3")
SET(WRAPPER_LIBRARY_LINK_LIBRARIES ITKCommon ITKStatistics)

WRAPPER_LIBRARY_CREATE_WRAP_FILES()
WRAPPER_LIBRARY_CREATE_LIBRARY()

