#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8.12)
cmake_policy(VERSION 2.8.12)

crag_library(CryptoKL
  KLKeyGeneraion
)

target_link_libraries(CryptoKL
  PUBLIC BraidGroup
)
