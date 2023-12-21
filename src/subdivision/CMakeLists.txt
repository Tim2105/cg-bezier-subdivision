file(GLOB_RECURSE SRCS ./*.cpp)
file(GLOB_RECURSE HDRS ./*.h)

add_executable(subdivision ${SRCS} ${HDRS})
target_link_libraries(subdivision pmp_vis imgui glfw glew)
