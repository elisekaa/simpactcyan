cmake_minimum_required(VERSION 2.8)
if (CMAKE_VERSION VERSION_LESS "2.8.12")
	message("Version is lower than 2.8.12, don't need to set CMP0022")
else()
	cmake_policy(SET CMP0022 OLD) # Don't show warning
endif()

project(simpact-cyan)
set(VERSION 1.0.0)
set(CMAKE_CXX_STANDARD 14)
include(${PROJECT_SOURCE_DIR}/cmake/SimpactMacros.cmake)

simpact_setup()

# This contains the main simpact program
add_subdirectory(src)

# Install everything from the data dir in e.g. /usr/local/share/simpact-cyan/
install(DIRECTORY ${PROJECT_SOURCE_DIR}/data/ DESTINATION ${CMAKE_INSTALL_PREFIX}/share/simpact-cyan/)
install(FILES ${PROJECT_SOURCE_DIR}/python/pysimpactcyan.py DESTINATION ${CMAKE_INSTALL_PREFIX}/share/simpact-cyan/python/)
