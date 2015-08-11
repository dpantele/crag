#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

crag_library(DehnDiagram)
SET_TARGET_PROPERTIES(DehnDiagram PROPERTIES LINKER_LANGUAGE CXX)

crag_test(dehn_diagram
  DehnDiagram
)
