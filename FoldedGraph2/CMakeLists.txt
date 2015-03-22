#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

crag_library(CompressedWord
  compressed_word
)

crag_library(DisjointSet)

crag_library(FoldedGraph2
  folded_graph
)

target_link_libraries(FoldedGraph2 CompressedWord)

crag_test(folded_graph 
  #AccHarvest
  FoldedGraph2
)

crag_test(cword CompressedWord)

crag_test(disjoint_set DisjointSet)

crag_test(xy_letter CompressedWord)

if(UNIX)
  target_compile_options(CompressedWord PUBLIC -Werror)
  target_compile_options(DisjointSet PUBLIC -Werror)
endif()
