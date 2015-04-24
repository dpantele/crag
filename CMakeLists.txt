#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 3.2)
cmake_policy(VERSION 3.2)
IF(UNIX)
  #add -std=c++11 for all subprojects to make ABI fully compatible
  #also enables eclipse c++11 features
  #dirty hack, but the only one I can find
  set(CMAKE_CXX_COMPILER_ARG1 "-std=c++11")
  add_definitions(-std=c++11)
ENDIF(UNIX)

# Name of the project
project(crag)

enable_testing()

if(BUILD_SHARED_LIBS)
  message(FATAL_ERROR "Shared build is not supported")
endif()

IF(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W3 -EHsc -nologo")
  foreach (flag_var
    CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    if(NOT BUILD_SHARED_LIBS)
      #force using static runtime
      string(REPLACE "/MD" "-MT" ${flag_var} "${${flag_var}}")
    endif()
  endforeach()
ENDIF(MSVC)

if(UNIX)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
endif(UNIX)

set_property(GLOBAL PROPERTY USE_FOLDERS ON) #turn on folders for VS

include("cmake/common.cmake")

file(GENERATE
  OUTPUT "${CMAKE_BINARY_DIR}/dumb.cpp"
  CONTENT "namespace {void dumb(){(void) dumb;}}"
)

add_library(dumb OBJECT "${CMAKE_BINARY_DIR}/dumb.cpp")

# Add subdirectories (All new directories should go here)
# add_subdirectory(3rdParty/ezOptionParser)
add_subdirectory(3rdParty/boost)
add_subdirectory(3rdParty/gtest)
add_subdirectory(3rdParty/gmp)
add_subdirectory(Alphabet)
add_subdirectory(BraidGroup)
add_subdirectory(CryptoAAG)
add_subdirectory(CryptoAE)
add_subdirectory(CryptoKL)
add_subdirectory(CryptoShftConj)
add_subdirectory(CryptoTripleDecomposition)
add_subdirectory(Elt)
add_subdirectory(Equation)
add_subdirectory(Examples)
add_subdirectory(FoldedGraph2)
add_subdirectory(FoldedGraphAccHarvest)
add_subdirectory(FreeGroup)
add_subdirectory(general)
add_subdirectory(Graph)
add_subdirectory(Graphics)
add_subdirectory(Group)
add_subdirectory(HigmanGroup)
add_subdirectory(Maps)
add_subdirectory(ranlib)
add_subdirectory(StringSimilarity)
add_subdirectory(SbgpFG)
add_subdirectory(TheGrigorchukGroup)
