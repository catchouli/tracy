project(tracy VERSION 0.1 LANGUAGES CXX)

add_executable(tracy
  private/main.cpp
)

enableCXX17(tracy)

target_link_libraries(tracy stb SDL2 glad glm SDL2 tiny_obj_loader simplefilewatcher imgui)

target_include_directories(tracy PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/public>
    $<INSTALL_INTERFACE:public>
    PRIVATE private)