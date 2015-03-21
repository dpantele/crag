#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

crag_library(Maps
  Map
  WhiteheadAutoSet
)

target_link_libraries(Maps
  PUBLIC Elt
  PRIVATE ranlib
)

#crag_main(test_rand Graph)
