cmake_minimum_required(VERSION 2.6)
add_executable(example example.cpp)
target_link_libraries(example openbabel)
target_link_libraries(example ~/Downloads/openbabel-2.3.2/openbabel-install/lib/libopenbabel.so)
include_directories(~/Downloads/openbabel-2.3.2/openbabel-install/include/openbabel-2.0)