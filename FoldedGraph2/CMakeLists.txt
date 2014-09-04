#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include("../cmake/common.cmake")

crag_library(FoldedGraph2
  folded_graph
)

crag_test(folded_graph FoldedGraph2)
crag_main(acc FoldedGraph2)