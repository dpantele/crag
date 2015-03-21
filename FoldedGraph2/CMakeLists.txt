#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

crag_library(FoldedGraph2
  folded_graph
  compressed_word
)

crag_test(folded_graph 
  AccHarvest
  FoldedGraph2
)

crag_test(cword FoldedGraph2)
