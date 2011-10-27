add_executable(tutorial tutorial.cpp)
target_link_libraries(tutorial ${GMP_LIBRARIES} ${GMPXX_LIBRARIES})