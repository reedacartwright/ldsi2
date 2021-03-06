CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)
PROJECT(ldsi)

SET(CMAKE_BUILD_TYPE_INIT "RelWithDebInfo")
SET(CMAKE_VERBOSE_MAKEFILE OFF CACHE BOOL "Turn on Verbose Makefiles" FORCE)
SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Modules")
SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)

SET(VERSION_MAJOR "2")
SET(VERSION_MINOR "0")
SET(VERSION_PATCH "0")
SET(VERSION_BRANCH "CURRENT")

INCLUDE(CheckIncludeFileCXX)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckLibraryExists)

FIND_PACKAGE(Git)
IF(GIT_FOUND)
  EXECUTE_PROCESS(
    COMMAND ${GIT_EXECUTABLE} describe --tags --dirty
      RESULT_VARIABLE GIT_RESULT
	  OUTPUT_VARIABLE PKG_VERSION
	  OUTPUT_STRIP_TRAILING_WHITESPACE
	  ERROR_QUIET
  )
  IF(GIT_RESULT EQUAL 0)
    SET(VERSION_REV "${PKG_VERSION}")
  ENDIF()
ENDIF()

SET(LDSI_VERSION "${VERSION_MAJOR}")
IF(${VERSION_MINOR})
  SET(LDSI_VERSION "${LDSI_VERSION}.${VERSION_MINOR}")
ENDIF(${VERSION_MINOR})
IF(${VERSION_PATCH})
  SET(LDSI_VERSION "${LDSI_VERSION}.${VERSION_PATCH}")
ENDIF(${VERSION_PATCH})
IF(DEFINED VERSION_BRANCH)
	SET(LDSI_VERSION "${LDSI_VERSION}-${VERSION_BRANCH}")
ENDIF(DEFINED VERSION_BRANCH)
IF(NOT ${VERSION_BRANCH} MATCHES RELEASE)
	SET(LDSI_VERSION "${LDSI_VERSION}-r${VERSION_REV}")
ENDIF(NOT ${VERSION_BRANCH} MATCHES RELEASE)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "LDSI")
SET(CPACK_PACKAGE_VENDOR "Reed A. Cartwright, PhD")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/readme.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/copying.txt")
SET(CPACK_PACKAGE_VERSION "${LDSI_VERSION}")
SET(CPACK_PACKAGE_VERSION_MAJOR "${VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${VERSION_PATCH}")

#Bug in cpack
SET(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME})
IF(${CPACK_SYSTEM_NAME} MATCHES Windows)
  IF(CMAKE_CL_64)
    SET(CPACK_SYSTEM_NAME win64)
  ELSE(CMAKE_CL_64)
    SET(CPACK_SYSTEM_NAME win32)
  ENDIF(CMAKE_CL_64)
ENDIF(${CPACK_SYSTEM_NAME} MATCHES Windows)

SET(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_SYSTEM_NAME}")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-Source")

SET(CPACK_PACKAGE_INSTALL_DIRECTORY "LDSI")

SET(CPACK_SOURCE_IGNORE_FILES
  "/CVS/"
  "/\\\\.svn/"
  "\\\\.swp$"
  "\\\\.#"
  "/#"
  ".*~$"
  "ldsi[-]"
  "/CMakeFiles/"
  "CMakeCache.txt"
  "CPack.*Config.cmake"
  "cmake_install.cmake"
  "_CPack_Packages"
)

IF(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
  #SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\\\InstallIcon.bmp")
  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\ldsi.exe")
  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} LDSI")
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\scit.us\\\\")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\scit.us\\\\")
  SET(CPACK_NSIS_CONTACT "reed@scit.us")
  SET(CPACK_NSIS_MODIFY_PATH ON)
  SET(CPACK_STRIP_FILES "bin/ldsi.exe")
ELSE(WIN32 AND NOT UNIX)
  SET(CPACK_GENERATOR "TGZ")
  SET(CPACK_SOURCE_GENERATOR "TGZ") 
  SET(CPACK_STRIP_FILES "bin/ldsi")
ENDIF(WIN32 AND NOT UNIX)

SET(CPACK_PACKAGE_EXECUTABLES "ldsi" "LDSI")
INCLUDE(CPack)

ADD_SUBDIRECTORY(src)

#INSTALL(FILES anakl RUNTIME DESTINATION bin)

#INCLUDE(InstallRequiredSystemLibraries)
 
