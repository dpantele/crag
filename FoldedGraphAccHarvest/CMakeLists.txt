#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

crag_library(AccHarvest
  acc
)

target_link_libraries(AccHarvest PUBLIC FoldedGraph2)
crag_test(acc FoldedGraph2 AccHarvest)

OPTION(FOLDED_GRAPH_ACC_TRY_SPLIT "If we are processing pair (u, v) = (as, bs), where |s| > |v| / 2, simply replace it with (aB, bs)" Off)

crag_main(main AccHarvest)
find_package(Threads REQUIRED)
target_link_libraries(FoldedGraphAccHarvest_main_main PRIVATE "${CMAKE_THREAD_LIBS_INIT}")
if(FOLDED_GRAPH_ACC_TRY_SPLIT)
  target_compile_definitions(FoldedGraphAccHarvest_main_main PRIVATE CRAG_FOLDED_GRAPH_ACC_TRY_SPLIT)
endif()

crag_main(auto_check AccHarvest)
crag_main(compare_generated_words AccHarvest)
crag_main(conj_classes AccHarvest)
crag_main(count_length AccHarvest)
crag_main(find_canonical AccHarvest)
crag_main(pairs_complexity AccHarvest)
crag_main(pick_random_pairs AccHarvest)
