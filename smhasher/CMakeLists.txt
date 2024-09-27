add_library(smhasher MurmurHash3.cpp)
target_include_directories(smhasher PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})