#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

include("../cmake/common.cmake")

crag_library(AccHarvest
  acc
)

target_link_libraries(AccHarvest PUBLIC FoldedGraph2)
crag_test(acc FoldedGraph2 AccHarvest)

crag_main(main AccHarvest)
find_package(Threads REQUIRED)
target_link_libraries(FoldedGraphAccHarvest_main_main PRIVATE "${CMAKE_THREAD_LIBS_INIT}")

crag_main(auto_check AccHarvest)
crag_main(compare_generated_words AccHarvest)
crag_main(count_length AccHarvest)
crag_main(pick_random_pairs AccHarvest)
