#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 3.1.0)
cmake_policy(VERSION 3.1.0)

crag_library(FoldedGraph2
  folded_graph
  compressed_word
)

target_compile_features(FoldedGraph2 PUBLIC cxx_generalized_initializers)

crag_test(folded_graph 
  AccHarvest
  FoldedGraph2
)

crag_test(cword FoldedGraph2)
