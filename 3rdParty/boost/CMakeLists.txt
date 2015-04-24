cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include(ExternalProject)

ExternalProject_Add(3rdPartyBoostBuild
  URL http://sourceforge.net/projects/boost/files/boost/1.58.0/boost_1_58_0.tar.bz2/download
  URL_HASH MD5=b8839650e61e9c1c0a89f371dd475546
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ${CMAKE_COMMAND} -E rename "<SOURCE_DIR>/boost" "<INSTALL_DIR>/include/boost"
)

ExternalProject_Get_Property(3rdPartyBoostBuild INSTALL_DIR)

file(MAKE_DIRECTORY "${INSTALL_DIR}/include")

add_library(3rdPartyBoost STATIC $<TARGET_OBJECTS:dumb>)
add_dependencies(3rdPartyBoost 3rdPartyBoostBuild)
add_library(3rdPartyBoostNoThreads STATIC $<TARGET_OBJECTS:dumb>)
add_dependencies(3rdPartyBoostNoThreads 3rdPartyBoostBuild)

set_target_properties(3rdPartyBoost 3rdPartyBoostNoThreads PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_DIR}/include"
  LINKER_LANGUAGE CXX
)

target_compile_definitions(3rdPartyBoostNoThreads INTERFACE BOOST_DISABLE_THREADS)

crag_test(test_boost_optional)
target_link_libraries(boost_test_test_boost_optional PRIVATE 3rdPartyBoost)