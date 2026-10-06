cmake_minimum_required(VERSION 3.16)

project(ImagineModelsConsumer LANGUAGES CXX)

find_package(ImagineModels 0.2 REQUIRED)

add_executable(consumer main.cc)
target_link_libraries(consumer PRIVATE ImagineModels::ImagineModels)
