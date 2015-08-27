#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 3.1.0)
cmake_policy(VERSION 3.1.0)

include("../cmake/common.cmake")

crag_library(AccHarvest
  acc
)

target_compile_features(AccHarvest PRIVATE cxx_range_for)

target_link_libraries(AccHarvest PUBLIC FoldedGraph2)
crag_test(acc FoldedGraph2 AccHarvest)

crag_main(auto_check AccHarvest)
crag_main(chains AccHarvest)
crag_main(main AccHarvest)

