#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

crag_library(boost_pool_gmp_allocator
  gmp_boost_pool_allocator
)

target_link_libraries(boost_pool_gmp_allocator
  PRIVATE 3rdPartyBoostNoThreads
  PRIVATE gmp
)

crag_library(crag_general
  Permutation
  ConfigFile
  BalancedTree
  VectorEnumerator
)

target_link_libraries(crag_general PRIVATE ranlib)

crag_main(mask crag_general ranlib)
